#pragma once

#include <string>
#include <vector>



/**
 *  Helper class for ultra-fast bootstrap
 *  Stores the best bs tree, its likelihood,
 *  and the bootstrap sampling
 */ 
class Bootstrap {
public:
  /**
   *  Elements is the number of elements in the distribution to sample
   *  for the local parallel core (the boostrapping will be done globally
   *  over all cores)
   */
  Bootstrap(unsigned int elements);

  /**
   * Evaluate the likelihood for this bootstrap from the per-sample likelihoods
   */
  double evaluate(const std::vector<double> &likelihoods) const;

private:
  std::vector<unsigned int> indices;
};



class RootBoot {
public:
  /**
   *  elements: see the class Bootstrap
   */
  RootBoot(unsigned int elements);

  /**
   *  Evaluate the bootstraped likelihood, and keep id as the best id if its
   *  value the best bootstraped likelihood encountered so far
   */
  void testRoot(const std::vector<double> &values, unsigned int id);
  
  /**
   *  Reset the best likelihood and best ID
   */
  void reset();

  /**
   *  Get the best ID encountered so far
   */
  unsigned int getBestID() const {return bestId;}
private: 
  Bootstrap bootstrap;
  unsigned int bestId;
  double bestLL;

};

class PerBranchBoot {
public:
  /**
   *  elements: see the class Bootstrap
   *  branches: number of branches
   */
  PerBranchBoot(unsigned int elements, unsigned int branches);

  /**
   *  Evaluate the bootstraped likelihood, and for each branch in branches: update
   *  the best likelihood so far and set isOk to isReferenceTree for this branch if
   *  the likelihood was updated
   */
  void test(const std::vector<double> &values, 
      const std::vector<unsigned int> &branches,
      bool isReferenceTree);

  /**
   *  Reset the best likelihoods and isOk values
   */
  void reset();

  /**
   *  Return true if this branch is ok (corresponds to the reference tree) 
   *  for this bootstrap
   */ 
  bool isOk(unsigned int branch) const {return _ok[branch];}
private:
  Bootstrap _bootstrap;
  std::vector<double> _bestLLs;
  std::vector<bool> _ok;

};


class PerBranchKH {
public:
  PerBranchKH(unsigned int elements, unsigned int branches, unsigned int bootstraps);

  void test(const std::vector<double> &values, 
      const std::vector<unsigned int> &branches);

  void newMLTree(const std::vector<double> &values);
  void newML(const std::vector<double> &values);

  unsigned int getSupport(unsigned int branch) const {return _oks[branch];}
private:
  std::vector<Bootstrap> _bootstraps;
  double _refLL;
  std::vector<double> _perBootstrapRefLL;
  std::vector<unsigned int> _oks;
};

