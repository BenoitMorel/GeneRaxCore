#include "SpeciesSearchCommon.hpp"

#include <trees/SpeciesTree.hpp>
#include <trees/PLLRootedTree.hpp>


static std::string getSubtreeNewick(const corax_rnode_t *subtree)
{
  if (!subtree->left) {
    return std::string(subtree->label);
  }
  std::string res("(");
  std::string id1 = getSubtreeNewick(subtree->left);
  std::string id2 = getSubtreeNewick(subtree->right);
  if  (id1 > id2) {
    std::swap(id1, id2);
  }
  return std::string("(") + id1 + "," + id2 + ")";
}
  
bool RootLikelihoods::hasSubtreeId(const corax_rnode_t *subtree) const {
  auto newickID = getSubtreeNewick(subtree);
  return newickToId.find(newickID) != newickToId.end();
}

unsigned int RootLikelihoods::getSubtreeId(const corax_rnode_t *subtree) const
{
  auto newickID = getSubtreeNewick(subtree);
  auto it = newickToId.find(newickID);
  return it->second;
}

unsigned int RootLikelihoods::getRootId(const corax_rnode_t *root)
{
  auto newickID1 = getSubtreeNewick(root->left);
  auto newickID2 = getSubtreeNewick(root->right);
  auto it1 = newickToId.find(newickID1);
  auto it2 = newickToId.find(newickID2);
  if (it1 == newickToId.end() && it2 == newickToId.end()) {
    unsigned int id = newickToId.size() / 2;
    newickToId.insert({newickID1, id});
    newickToId.insert({newickID2, id});
    return id;
  }
  assert(it1->second == it2->second);
  return it1->second;
}

void RootLikelihoods::saveRootLikelihood(corax_rnode_t *root, double ll) 
{
  auto id = getRootId(root);
  _idToLL[id] = ll;
}
  
void RootLikelihoods::savePerFamilyLikelihoods(corax_rnode_t *root, 
      const PerFamLL &likelihoods)
{
  auto id = getRootId(root);
  for (auto &bs: _bootstraps) {
    bs.testRoot(likelihoods, id);  
  }
}

void RootLikelihoods::fillTree(PLLRootedTree &tree)
{
  std::vector<double> nodeIdToLL(tree.getNodeNumber(), 0.0);
  double bestLL = -std::numeric_limits<double>::infinity();
  for (auto node: tree.getNodes()) {
    if (!hasSubtreeId(node)) {
      continue;
    }
    // we have a likelihood value
    auto id = getSubtreeId(node);
    auto value = _idToLL[id];
    nodeIdToLL[node->node_index] = value;
    bestLL = std::max<double>(value, bestLL);
  }
  for (auto node: tree.getNodes()) {
    bool hasValue = nodeIdToLL[node->node_index] != 0.0;
    std::string label;
    if (hasValue) {
      double value = nodeIdToLL[node->node_index] - bestLL;
      if(!node->left && node->label) {
        label = std::string(node->label);
        label += "_";
      }
      label += std::to_string(value);
      tree.setLabel(node->node_index, label);
    }
  }
}

void RootLikelihoods::fillTreeBootstraps(PLLRootedTree &tree)
{
  std::unordered_map<unsigned int, unsigned int> idToSupport;
  for (auto it: _idToLL) {
    idToSupport.insert({it.first, 0});
  }
  for (const auto &bs: _bootstraps) {
    idToSupport[bs.getBestID()] ++;
  }
  auto postOrderNodes = tree.getPostOrderNodes();
  for (auto it = postOrderNodes.rbegin(); it != postOrderNodes.rend(); ++it) {
    
    auto node = *it;
    if (!hasSubtreeId(node)) {
      continue;
    }
    auto id = getSubtreeId(node);
    if (idToSupport.find(id) == idToSupport.end()) {
      continue;
    }
    std::string label;
    auto value = idToSupport[id];
    if(!node->left && node->label) {

      label = std::string(node->label);
      free(node->label);
      node->label = nullptr;
      label += "_";
    }
    label += std::to_string(value);
    node->label = (char*)malloc(label.size() + 1);
    memcpy(node->label, label.c_str(), label.size());
    node->label[label.size()] = 0;
  }
}

void SpeciesSearchState::betterTreeCallback(double ll, PerFamLL &perFamLL)
{
  bool masterRankOnly = true;
  speciesTree.saveToFile(pathToBestSpeciesTree, masterRankOnly);
  bestLL = ll;
  khBoots.newMLTree(perFamLL);
  for (auto listener: _listeners) {
    listener->betterTreeCallback();
  }
}

void SpeciesSearchState::betterLikelihoodCallback(double ll, PerFamLL &perFamLL)
{
  bestLL = ll;
  khBoots.newML(perFamLL);
}

bool SpeciesSearchCommon::testSPR(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
    SpeciesSearchState &searchState,
    unsigned int prune,
    unsigned int regraft
    )
{
  evaluation.pushRollback();
  // Apply the move
  std::vector<unsigned int> affectedBranches;
  SpeciesTreeOperator::getAffectedBranches(speciesTree, prune, regraft,
      affectedBranches);
  auto rollback = SpeciesTreeOperator::applySPRMove(speciesTree, prune, regraft);
  bool runExactTest = true;
  double approxLL = 0.0;
  if (evaluation.providesFastLikelihoodImpl()) {
    // first test with approximative likelihood
    approxLL = evaluation.computeLikelihoodFast();
    if (searchState.averageApproxError.isSignificant()) {
      //  Decide whether we can already
      // discard the move
      auto epsilon = 2.0 * searchState.averageApproxError.getAverage();
      runExactTest &= (approxLL + epsilon > searchState.bestLL );
    } 
  }
  if (runExactTest) {
    // we test the move with exact likelihood
    PerFamLL perFamLL;
    auto testedTreeLL = evaluation.computeLikelihood(&perFamLL);
    for (auto &bs: searchState.sprBoots) {
      bs.test(perFamLL, affectedBranches, false); 
    }
    if (evaluation.providesFastLikelihoodImpl()) {
      searchState.averageApproxError.addValue(testedTreeLL - approxLL);
    }
    if (testedTreeLL > searchState.bestLL + 0.00000001) {
      searchState.betterTreeCallback(testedTreeLL, perFamLL);
      // Better tree found! Do not rollback, and return
      return true;
    } else {
      searchState.khBoots.test(perFamLL, affectedBranches);
    }
  }
  // the tree is not better, rollback the move
  SpeciesTreeOperator::reverseSPRMove(speciesTree, prune, rollback);
  evaluation.popAndApplyRollback();
  return false;
}

bool SpeciesSearchCommon::veryLocalSearch(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
    SpeciesSearchState &searchState,
    unsigned int spid)
{

  const unsigned int radius = 2;
  std::vector<unsigned int> prunes;
  prunes.push_back(spid);
  for (auto prune: prunes) {
    std::vector<unsigned int> regrafts;
    SpeciesTreeOperator::getPossibleRegrafts(speciesTree, 
        prune, 
        radius, 
        regrafts);
    for (auto regraft: regrafts) {
      if (testSPR(speciesTree, evaluation, searchState, prune, regraft)) {
        Logger::timed << "\tfound better* (LL=" 
            << searchState.bestLL << ", hash=" << 
            speciesTree.getHash() << ")" << std::endl;
        veryLocalSearch(speciesTree, evaluation, searchState,
            prune);
        return true;
      }
    }
  }
  return false;
}

void SpeciesSearchState::saveSpeciesTreeBP(const std::string &outputFile)
{
  PLLRootedTree tree(speciesTree.getTree().getNewickString(), false);
  auto m = tree.getNodeIndexMapping(speciesTree.getTree());  
  for (auto node: tree.getNodes()) {
    if (!node->left) {
      continue;
    }
    unsigned int ok = 0;
    for (const auto &bs: sprBoots) {
      if (bs.isOk(m[node->node_index])) {
        ok++;
      }
    }
    auto label = std::to_string(ok);
    tree.setLabel(node->node_index, label);
  }
  tree.save(outputFile);

}

void SpeciesSearchState::saveSpeciesTreeKH(const std::string &outputFile)
{
  PLLRootedTree tree(speciesTree.getTree().getNewickString(), false);
  auto m = tree.getNodeIndexMapping(speciesTree.getTree());  
  for (auto node: tree.getNodes()) {
    if (!node->left) {
      continue;
    }
    auto ok = khBoots.getSupport(m[node->node_index]);
    auto label = std::to_string(ok);
    tree.setLabel(node->node_index, label);
  }
  tree.save(outputFile);

}

