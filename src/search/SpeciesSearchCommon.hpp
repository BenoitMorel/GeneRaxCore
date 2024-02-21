#pragma once

#include <vector>

#include <likelihoods/ReconciliationEvaluation.hpp>
#include <util/types.hpp>
#include <util/Scenario.hpp>
#include <maths/AverageStream.hpp>
#include <search/UFBoot.hpp>
#include <trees/SpeciesTree.hpp>

class PerCorePotentialTransfers;
using TreePerFamLL = std::pair<std::string, PerFamLL>;
using TreePerFamLLVec = std::vector<TreePerFamLL>;


/**
 *  Store results (likelihoods, bootstrap info) for each
 *  root candidate that has been evaluated
 */
class RootLikelihoods {
public:

  RootLikelihoods(unsigned int localFamilies)
  {
    unsigned int samples = 1000;
    for (unsigned int i = 0; i < samples; ++i) {
      _bootstraps.push_back(RootBoot(localFamilies));
    }
  }
                

  /*
   *  Reset all results
   */
  void reset() {
    _idToLL.clear();
    newickToId.clear();
    for (auto &bs: _bootstraps) {
      bs.reset();
    }
  }

  /**
   *  Add the likelihood of a root position
   */
  void saveRootLikelihood(corax_rnode_t *root, double ll);
  
  void savePerFamilyLikelihoods(corax_rnode_t *root, 
      const PerFamLL &likelihoods);

  /**
   *  Fill the labels of the tree with the corresponding root 
   *  likelihoods 
   */
  void fillTree(PLLRootedTree &tree);

  /**
   *  Fill the labels of the tree with the corresponding root 
   *  likelihoods 
   */
  void fillTreeBootstraps(PLLRootedTree &tree);

  /**
   * Return true if no value was stored
   */ 
  bool isEmpty() const {return !_idToLL.size();}

private:

  bool hasSubtreeId(const corax_rnode_t *subtree) const;
  unsigned int getSubtreeId(const corax_rnode_t *subtree) const;
  unsigned int getRootId(const corax_rnode_t *root);

  std::unordered_map<std::string, unsigned int> newickToId;
  std::unordered_map<unsigned int, double> _idToLL;
  std::vector<RootBoot> _bootstraps;
};

/**
 *  Interface for classes that are used by the search
 *  algorithms to evaluate the reconciliation likelihood and 
 *  to describe the reconciliation model.
 */
class SpeciesTreeLikelihoodEvaluatorInterface {
public:
  virtual ~SpeciesTreeLikelihoodEvaluatorInterface() {};
  /**
   * Computes and return the likelihood
   *
   * If perFam is set, fill perFamLL with the per-family log-likelihoods
   * from the current parralel core 
   */
  virtual double computeLikelihood(PerFamLL *perFamLL = nullptr) = 0;

  /**
   *  Fast but approximated version of the likelihood
   *  computation
   */
  virtual double computeLikelihoodFast() = 0;

  /**
   *  Return true if computeLikelihood and computeLikelihoodFast
   *  have different implementations
   */
  virtual bool providesFastLikelihoodImpl() const = 0;

  /**
   *  Return true if the model is dated (if the model depends
   *  on the speciation event order)
   */
  virtual bool isDated() const = 0;

  /**
   *  Optimize model rates, such as DTL rates.
   */
  virtual double optimizeModelRates(bool thorough = false) = 0;

  /**
   *  Save in a stack what needs to be saved in case 
   *  of the rollback of a species tree operation
   */
  virtual void pushRollback() = 0;

  /**
   * Pop the upper state and apply it after a rollback
   * of a species tree operation
   */
  virtual void popAndApplyRollback() = 0;
  
  virtual void getTransferInformation(SpeciesTree &speciesTree,
    TransferFrequencies &frequencies,
    PerSpeciesEvents &perSpeciesEvents,
    PerCorePotentialTransfers &potentialTransfers) = 0;
  
  
  /**
   *  Are we in prune species tree mode?
   */
  virtual bool pruneSpeciesTree() const = 0;

  /**
   *  Should be called when the species tree is updated
   */
  virtual void onSpeciesTreeChange(
      const std::unordered_set<corax_rnode_t *> *nodesToInvalidate) {
    (void)(nodesToInvalidate);
  }

  /**
   *  Should be called when the species tree dates (speciation orders)
   *  change
   */
  virtual void onSpeciesDatesChange() {};
};

/**
 *  Structure that describes the current state of the
 *  species tree search
 */
struct SpeciesSearchState {
public:
  SpeciesSearchState(SpeciesTree &speciesTree,
      const std::string &pathToBestSpeciesTree,
      unsigned int familyNumber): 
    speciesTree(speciesTree),
    pathToBestSpeciesTree(pathToBestSpeciesTree),
    farFromPlausible(true),
    khBoots(familyNumber, speciesTree.getTree().getNodeNumber(), 1000)
    {
      for (unsigned int i = 0; i < 1000; ++i) {
        sprBoots.push_back(PerBranchBoot(familyNumber,
              speciesTree.getTree().getNodeNumber()));
      }
    }
  
  /**
   *  Reference to the current species tree
   */
  SpeciesTree &speciesTree;
  
  /**
   *  The search algorithm will save the current species tree
   *  at this location after each likelihood improvement
   */
  std::string pathToBestSpeciesTree;

  /**
   *  The likelihood of the best species tree
   */
  double bestLL;
        
   /**
   *  Set to true when the tree is far from being plausible
   *  (in the early stage of the search after starting from
   *  a random tree for instance).
   *  
   *  When set to true, the search strategy can choose to
   *  apply less thorough optimizations and to favor speed
   *  over small tree improvments.
   *
   *  Both the search algorithm functions and their caller 
   *  can update this variable.
   */
  bool farFromPlausible;

  /**
   *  Stores the average difference between the average and
   *  the exact likelihood values. This average is updated
   *  in a streaming fashion everytime that both the average
   *  and exact likelihood are computed during the search
   *  (see for instance SpeciesSearchCommon::testSPR)
   *  It is then used to decide if the approximated likelihood
   *  score is good enough to try estimating the exact score.
   *  
   *  This is only relevant when 
   *  SpeciesTreeLikelihoodEvaluatorInterface::providesFastLikelihoodImpl
   *  is set to true
   */
  AverageStream averageApproxError;


  std::vector<PerBranchBoot> sprBoots;
  PerBranchKH khBoots;  

  /**
   *  To call when a better tree is found
   */
  void betterTreeCallback(double ll, PerFamLL &perFamLL);

  /**
   *  To call when the likelihood increases but the tree does not 
   *  change (e.g after rates optimization)
   */
  void betterLikelihoodCallback(double ll, PerFamLL &perFamLL);

  void saveSpeciesTreeKH(const std::string &outputFile);
  void saveSpeciesTreeBP(const std::string &outputFile);
  
  class Listener {
  public:
    /**
     *  If the object is attached to a SpeciesSearchState, 
     *  this callback will be called when a better tree is found
     */
    virtual void betterTreeCallback() = 0;
  };

  void addListener(Listener *listener) {_listeners.push_back(listener);}

private:
  std::vector<Listener *> _listeners;
};

class SpeciesSearchCommon {
public:
  /**
   *  Test a SPR move. 
   *  If it improves the likelihood, keep it and return true.
   *  Else, rollback it and return false.
   */
  static bool testSPR(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
    SpeciesSearchState &searchState,
    unsigned int prune,
    unsigned int regraft);

  /**
   *  Try SPR moves with a small radius around the species
   *  node with id spid. If a move improves the likelihood, 
   *  apply it, and recursively search further. Otherwise, 
   *  rollback the tested moves.
   *  Returns true if one better tree has been found
   */
  static bool veryLocalSearch(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
    SpeciesSearchState &searchState,
    unsigned int spid);
  
};

