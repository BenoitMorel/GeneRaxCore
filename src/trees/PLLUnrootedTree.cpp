#include "PLLUnrootedTree.hpp"

#include <IO/Logger.hpp>
#include <trees/PLLRootedTree.hpp>  
#include <stack>
#include <functional>
#include <sstream>
#include <deque>
#include <IO/LibpllException.hpp>
#include <corax/tree/utree_compare.h>
#include <corax/io/newick.hpp>
#include <fstream>

void defaultUnodePrinter(corax_unode_t *node, 
    std::stringstream &ss)
{
  if (node->label) {
    ss << node->label;
  }
  ss << ":" << node->length;
}


static void destroyNodeData(void *)
{
}

void utreeDestroy(corax_utree_t *utree) {
  if(!utree)
    return;
  corax_utree_destroy(utree, destroyNodeData);
}



static corax_utree_t *readNewickFromStr(const std::string &str) 
{
  corax_newick_parser_t parser(str);
  auto utree = parser.parse(true, true);  
  //auto utree = corax_utree_parse_newick_string_unroot(str.c_str());
  if (!utree) 
    throw LibpllException("Error while reading tree from std::string: ", str);
  return utree;
}


static corax_utree_t *readNewickFromFile(const std::string &str)
{
  std::ifstream is(str);
  if (!is) {
    throw LibpllException("Can't open file: ", str);
  }
  std::string newick;
  if (!std::getline(is, newick)) {
    throw LibpllException("Error while reading tree from file: ", str);
  }
  return readNewickFromStr(newick);
}

static corax_utree_t *buildUtree(const std::string &str, bool isFile)
{
  if (isFile) {
    return readNewickFromFile(str);
  } else {
    return readNewickFromStr(str);
  }
}

PLLUnrootedTree::PLLUnrootedTree(const std::string &str, bool isFile):
  _tree(buildUtree(str, isFile), utreeDestroy)
{
  setMissingBranchLengths();
}

PLLUnrootedTree::PLLUnrootedTree(PLLRootedTree &rootedTree):
  _tree(corax_rtree_unroot(rootedTree.getRawPtr()), utreeDestroy)
{
  corax_unode_t *root = 
    _tree->nodes[_tree->tip_count + _tree->inner_count - 1];
  corax_utree_reset_template_indices(root, _tree->tip_count);
}
  
std::unique_ptr<PLLUnrootedTree> PLLUnrootedTree::buildFromStrOrFile(const std::string &strOrFile)
{
  std::unique_ptr<PLLUnrootedTree> res;
  try {
    res = std::make_unique<PLLUnrootedTree>(strOrFile, true);
  } catch (...) {
    res = std::make_unique<PLLUnrootedTree>(strOrFile, false);
  }
  return res;
}
  
  
std::unique_ptr<PLLUnrootedTree> PLLUnrootedTree::buildWithOutgroupFromRootedNewick(const std::string &newickRootedStr, const std::string &outgroup)
{
  std::string left = newickRootedStr;
  auto begin = left.find_first_of('(') + 1;
  auto end = left.find_last_of(')') - 1;
  left = left.substr(begin, end); 
  std::string newick = std::string("(") + left + "," + outgroup + ");";
  return std::make_unique<PLLUnrootedTree>(newick, false);
}
 
std::unique_ptr<PLLUnrootedTree> 
PLLUnrootedTree::buildWithOutgroupFromRooted(
    const PLLRootedTree &rootedTree, const std::string &outgroup)
{
  return buildWithOutgroupFromRootedNewick(rootedTree.getNewickString(),
      outgroup);
}

char *consensus_serialize(const corax_unode_t *node)
{
  
  std::string newick;
  assert(node->next != node);
  if (node->next) {
    
    auto data = static_cast<corax_consensus_data_t *>(node->data);
    if (data) {
      newick += std::to_string(data->support);
    }
  } else if (node->label) {
    newick += node->label;
  }
  // we don't have the branch lengths
  newick += ":0.1";
  char *res = static_cast<char*>(malloc(sizeof(char) * (newick.size() + 1)));
  res[newick.size()] = '\0';
  memcpy(res, newick.c_str(), newick.size());
  return res;
}

static void setNodeLabel(corax_unode_t *node, const std::string &label)
{
  free(node->label);
  node->label = static_cast<char*>(malloc(sizeof(char) * (label.size() + 1)));
  std::strcpy(node->label, label.c_str());
}

void PLLUnrootedTree::setLabel(unsigned int nodeIndex, const std::string &label)
{
  setNodeLabel(getNode(nodeIndex), label);
}

std::string PLLUnrootedTree::buildConsensusTree(
    std::vector<std::string> &strOrFiles, 
      double threshold)
{
  std::vector<std::shared_ptr<PLLUnrootedTree> >trees;
  for (const auto &str: strOrFiles) {
    trees.push_back(buildFromStrOrFile(str));
  }
  return buildConsensusTree(trees, threshold);
}

std::string PLLUnrootedTree::buildConsensusTree(
    const std::vector<std::shared_ptr<PLLUnrootedTree> > &trees,
      double threshold)
{
  assert(threshold >= 0.0 && threshold <= 1.0);
  std::vector<const corax_utree_t*> treePointers;
  std::vector<double> weights;
  for (auto &tree: trees) {
    treePointers.push_back(tree->getRawPtr());
    corax_utree_consistency_set(trees[0]->getRawPtr(), tree->getRawPtr());
  }
  auto weight = 1.0 / static_cast<double>(treePointers.size());
  weights = std::vector<double>(treePointers.size(), weight);
  assert(treePointers.size());
  assert(weights.size());
  auto consensus = corax_utree_weight_consensus(&(treePointers[0]),
        &(weights[0]),
        threshold,
        treePointers.size());
  assert(consensus);
  assert(consensus->tree);
  auto node = consensus->tree;
  auto newick = corax_utree_export_newick(node, consensus_serialize);
  return newick;
}


PLLUnrootedTree::PLLUnrootedTree(const std::vector<const char*> &labels,
    unsigned int seed):
  _tree(corax_utree_random_create(static_cast<unsigned int>(labels.size()), &labels[0], seed), utreeDestroy)
{

}

void PLLUnrootedTree::save(const std::string &fileName)
{
  std::ofstream os(fileName, std::ofstream::out);
  char *newick = corax_utree_export_newick_rooted(getRawPtr()->nodes[0], 0);
  os << newick;
  os.close();
  free(newick);
}

void PLLUnrootedTree::setMissingBranchLengths(double minBL)
{
  for (auto node: getLeaves()) {
    if (0.0 >= node->length) {
      node->length = minBL;
    } 
  }
  for (unsigned int i = _tree->tip_count; i < _tree->tip_count + _tree->inner_count; ++i) {
    if (0.0 >= _tree->nodes[i]->length)
      _tree->nodes[i]->length = minBL;
    if (0.0 >= _tree->nodes[i]->next->length)
      _tree->nodes[i]->next->length = minBL;
    if (0.0 >= _tree->nodes[i]->next->next->length)
      _tree->nodes[i]->next->next->length = minBL;
  }  
}
  
CArrayRange<corax_unode_t*> PLLUnrootedTree::getLeaves() const
{
  return CArrayRange<corax_unode_t*>(_tree->nodes, getLeafNumber());
}

CArrayRange<corax_unode_t*> PLLUnrootedTree::getNodes() const
{
  return CArrayRange<corax_unode_t*>(_tree->nodes, getNodeNumber());
}

CArrayRange<corax_unode_t*> PLLUnrootedTree::getInnerNodes() const
{
  return CArrayRange<corax_unode_t*>(_tree->nodes + getLeafNumber(), getInnerNodeNumber());
}


unsigned int PLLUnrootedTree::getNodeNumber() const
{
  return getLeafNumber() + getInnerNodeNumber();
}

unsigned int PLLUnrootedTree::getDirectedNodeNumber() const
{
  return _tree->edge_count * 2;
}

unsigned int PLLUnrootedTree::getLeafNumber() const
{
  return _tree->tip_count;
}

unsigned int PLLUnrootedTree::getInnerNodeNumber() const
{
  return _tree->inner_count;
}
  
corax_unode_t *PLLUnrootedTree::getNode(unsigned int node_index) const
{
  return _tree->nodes[node_index];
}

corax_unode_t *PLLUnrootedTree::getAnyInnerNode() const
{
  return getNode(getLeafNumber());
}
  
std::unordered_set<std::string> PLLUnrootedTree::getLeafLabels()
{
  std::unordered_set<std::string> res;
  for (auto leaf: getLeaves()) {
    if (leaf->label) {
      res.insert(std::string(leaf->label));
    }
  }
  return res;
}


static bool isBranchIn(corax_unode_t *b, 
    const std::unordered_set<corax_unode_t *> &branches)
{
  return branches.find(b) != branches.end() 
    || branches.find(b->back) != branches.end();
}

std::unordered_set<corax_unode_t *> PLLUnrootedTree::getBranches() const
{
  std::unordered_set<corax_unode_t *> branches;
  for (auto node: getNodes()) {
    if (!isBranchIn(node, branches)) {
      branches.insert(node);
    }
    if (node->next) {
      node = node->next;
      if (!isBranchIn(node, branches)) {
        branches.insert(node);
      }
      node = node->next;
      if (!isBranchIn(node, branches)) {
        branches.insert(node);
      }
    }
  }
  return branches;
}

struct ToSort {
  corax_unode_t *node;
  std::string label;
  unsigned int distanceToLeaf;
 
  ToSort() {}
  ToSort(corax_unode_t *inode, 
      const std::string &ilabel,
      unsigned int idistanceToLeaf): node(inode),
  label(ilabel),
  distanceToLeaf(idistanceToLeaf) {}
  ToSort(const ToSort &s): node(s.node),
    label(s.label),
    distanceToLeaf(s.distanceToLeaf) {}
};

static void fillBranchesRec(corax_unode_t * node,
    std::vector<ToSort> &toSortVector,
    ToSort &toSort)
{
  if (!node->next) {
    toSort = ToSort(node, std::string(node->label), 0);
    toSortVector.push_back(toSort);
    return;
  }
  ToSort toSort1;
  fillBranchesRec(node->next->back, toSortVector, toSort1);
  ToSort toSort2;
  fillBranchesRec(node->next->next->back, toSortVector, toSort2);
  if (toSort1.label > toSort2.label) {
    toSort = ToSort(node, toSort1.label, toSort1.distanceToLeaf + 1);
  } else {
    toSort = ToSort(node, toSort2.label, toSort2.distanceToLeaf + 1);
  }
  toSortVector.push_back(toSort);
}

struct less_than_key
{
  inline bool operator() (const ToSort& t1, const ToSort& t2)
  {
    if (t1.label == t2.label) {
      return t1.distanceToLeaf < t2.distanceToLeaf;
    } 
    return t1.label < t2.label;
  }
};

std::vector<corax_unode_t *> PLLUnrootedTree::getBranchesDeterministic() const
{
  // find a deterministic root
  corax_unode_t *root = getNode(0);
  std::string rootLabel(root->label);
  for (auto leaf: getLeaves()) {
    std::string label(leaf->label);
    if (rootLabel > label) {
      rootLabel = label;
      root = leaf;
    }
  }
  std::vector<ToSort> toSortVector;
  ToSort toSort;
  fillBranchesRec(root->back, toSortVector, toSort);
  assert(toSortVector.size() == getLeafNumber() * 2 - 3);
  std::sort(toSortVector.begin(), toSortVector.end(), less_than_key());
  std::vector<corax_unode_t *> res;
  for (const auto &toSort: toSortVector) {
    res.push_back(toSort.node);
  }
  return res;
}

static void fillPostOrder(corax_unode_t *node,
    std::vector<corax_unode_t*> &nodes,
    std::vector<char> &markedNodes)
{
  // we already traversed this node
  if (markedNodes[node->node_index]) {
    return;
  }
  // mark the node as traversed
  markedNodes[node->node_index] = true;
  // first process children
  if (node->next) {
    fillPostOrder(node->next->back, nodes, markedNodes);
    fillPostOrder(node->next->next->back, nodes, markedNodes);
  }
  nodes.push_back(node);
}
  
std::vector<corax_unode_t*> PLLUnrootedTree::getPostOrderNodesFrom(corax_unode_t *node) const
{
  std::vector<corax_unode_t*> nodes;
  std::vector<char> markedNodes(getDirectedNodeNumber(), false);
  fillPostOrder(node, nodes, markedNodes);
  return nodes;
}


std::vector<corax_unode_t*> PLLUnrootedTree::getPostOrderNodesRooted(corax_unode_t *virtualRoot) const
{
  std::vector<corax_unode_t*> nodes;
  std::vector<char> markedNodes(getDirectedNodeNumber(), false);
  fillPostOrder(virtualRoot, nodes, markedNodes);
  fillPostOrder(virtualRoot->back, nodes, markedNodes);
  return nodes;
}


std::vector<corax_unode_t*> PLLUnrootedTree::getPostOrderNodes(bool innerOnly) const
{
  std::vector<corax_unode_t*> nodes;
  std::vector<char> markedNodes(getDirectedNodeNumber(), false);
  if (innerOnly) {
    for (auto node: getLeaves()) {
      markedNodes[node->node_index] = true;
    }
  }
  // do the post order traversal from all possible virtual roots 
  for (auto node: getLeaves()) {
    fillPostOrder(node->back, nodes, markedNodes);
  }
  if (innerOnly) {
    assert(nodes.size() == getDirectedNodeNumber() - getLeafNumber());
  } else {
    assert(nodes.size() == getDirectedNodeNumber());
  }
  return nodes;
}
  
std::vector<corax_unode_t*> PLLUnrootedTree::getReverseDepthNodes() const
{
  std::deque<corax_unode_t *> q;
  std::vector<bool> marked(getDirectedNodeNumber(), false);
  std::vector<corax_unode_t *> nodes;
  for (auto leaf: getLeaves()) {
    marked[leaf->node_index] = true;
    q.push_back(leaf);
  }
  while (!q.empty()) {
    auto node = q.front();
    q.pop_front();
    nodes.push_back(node);
    auto back = node->back;
    if (!back->next) {
      continue;
    }
    auto left = back->next;;
    auto right = back->next->next;
    if (!marked[left->node_index]) {
      q.push_back(left);
      marked[left->node_index] = true;
    }
    if (!marked[right->node_index]) {
      q.push_back(right);
      marked[right->node_index] = true;
    }
  }
  assert(nodes.size() == getDirectedNodeNumber());
  return nodes; 
}

static void computePairwiseDistancesRec(corax_unode_t *currentNode, 
    double currentDistance,
    VectorDouble &distancesVector)
{
  currentDistance += currentNode->length;
  if (!currentNode->next) {
    // leaf
    distancesVector[currentNode->node_index] = currentDistance;
    return;
  }
  computePairwiseDistancesRec(currentNode->next->back, 
      currentDistance, 
      distancesVector);
  computePairwiseDistancesRec(currentNode->next->next->back, 
      currentDistance, 
      distancesVector);
} 

void PLLUnrootedTree::computePairwiseDistances(MatrixDouble &distances,
    bool leavesOnly)
{
  auto M = leavesOnly ? getLeafNumber() : getDirectedNodeNumber(); 
  auto N = getLeafNumber();
  auto MIter = leavesOnly ? getLeafNumber() : getNodeNumber();
  VectorDouble zeros(N, 0.0);
  distances = MatrixDouble(M, zeros);
  for (unsigned int i = 0; i < MIter; ++i) {
    auto node = getNode(i);
    auto &distancesVector = distances[node->node_index];
    computePairwiseDistancesRec(node->back, 0.0, distancesVector);
    if (node->next) {
      // compute distances to leaves in all three directions
      computePairwiseDistancesRec(node->next->back, 
          0.0, 
          distancesVector);
      computePairwiseDistancesRec(node->next->next->back, 
          0.0, 
          distancesVector);
      // also update the two other directed nodes
      distances[node->next->node_index] = distancesVector;
      distances[node->next->next->node_index] = distancesVector;
    }
  }
}


static void getCladeRec(corax_unode_t *node, 
    std::unordered_set<unsigned int> &clade)
{
  if (node->next) {
    getCladeRec(node->next->back, clade);
    getCladeRec(node->next->next->back, clade);
  } else {
    clade.insert(node->node_index);
  }
}

std::unordered_set<unsigned int> 
  PLLUnrootedTree::getClade(corax_unode_t *node)
{
  std::unordered_set<unsigned int> clade;
  getCladeRec(node, clade);
  return clade;
}


// look for *pv (or one of his nexts) under the oriented node u.
// if it is found, *pv is updated with the found next, and 
// the function returns true (and false otherwise)
static bool orientAux(corax_unode_t *u, 
    corax_unode_t **pv,
    std::stack<corax_unode_t *> &path)
{
  path.push(u);
  auto v = *pv;
  if (v == u) {
    return true;
  }
  if (v->next) {
    if (v->next == u) {
      *pv = v->next;
      return true;
    } else if (v->next->next == u) {
      *pv = v->next->next;
      return true;
    }
  }
  if (!u->next) { 
    // end of recursion, we did not find *v
    path.pop();
    return false;
  } else {
    if (orientAux(u->next->back, pv, path) || 
      orientAux(u->next->next->back, pv, path)) {
      return true;
    } else {
      path.pop();
      return false;
    }
  }
}

static void stackToVector(std::stack<corax_unode_t *> s,
    std::vector<corax_unode_t *> &v)
{
  v.clear();
  v.resize(s.size());
  for (int i = s.size() - 1; i >= 0; --i) {
    v[i] = s.top();
    s.pop();
  }
}

void PLLUnrootedTree::orientTowardEachOther(corax_unode_t **pu,
    corax_unode_t **pv,
    std::vector<corax_unode_t *> &branchesPath)
{
  assert((*pu) != (*pv));
  auto *u = *pu;
  std::stack<corax_unode_t *> path;
  if (orientAux(u->back, pv, path)) {
    stackToVector(path, branchesPath);    
    return;
  }
  assert(path.size() == 0);
  if (u->next) {
    if (orientAux(u->next->back, pv, path)) {
      *pu = u->next;
      stackToVector(path, branchesPath);    
      return;
    } else if (orientAux(u->next->next->back, pv, path)) {
      *pu = u->next->next;
      stackToVector(path, branchesPath);    
      return;
    }
  }
  assert(false);
}



std::vector<double> PLLUnrootedTree::getMADRelativeDeviations()
{
  MatrixDouble distances;
  computePairwiseDistances(distances, false);
  auto nodes = getPostOrderNodes();  
  std::vector<double> deviations(nodes.size());
  for (auto node: getPostOrderNodes()) {
    // we reuse notations from the MAD paper (I, J, i, j, rho, b, c)
    // I and J are the set of leaves of each side of the branch node
    // rho is the relative position of the rooting that minimizes 
    // the squared relative deviation
    auto I = getClade(node); 
    auto J = getClade(node->back);
    auto rho = 0.0;
    auto rhoDen = 0.0;
    auto i = node->node_index;
    auto Dij = node->length;
    for (auto b: I) {
      for (auto c: J) {
        auto invPowBC = pow(distances[b][c], -2.0);
        rho += (distances[b][c] - 2.0 * distances[i][b]) * invPowBC;
        rhoDen += invPowBC;
      }
    }
    rho = rho / (2.0 * Dij * rhoDen);
    rho = std::max(std::min(rho, 1.0), 0.0);
    auto deviation = 0.0;
    for (auto b: I) {
      for (auto c: J) {
        auto v = (2.0 * (distances[i][b] + rho * Dij) / distances[b][c]);
        v -= 1.0;
        deviation += pow(v, 2.0);
      }
    }
    deviations[node->node_index] = deviation;
  }
  return deviations;
}

static void printAux(corax_unode_t *node,
    std::stringstream &ss,
    UnodePrinter f)
{
  if (node->next) {
    ss << "(";
    auto temp = node->next;
    printAux(temp->back, ss, f);
    while (temp->next != node) {
      temp = temp->next;
      ss << ",";
      printAux(temp->back, ss, f);
    }
    ss << ")";
  }
  f(node, ss);
}
  
std::string PLLUnrootedTree::getSubtreeString(corax_unode_t *subtree, UnodePrinter f)
{
  if (!subtree) {
    return "null";
  }
  std::stringstream ss;
  printAux(subtree, ss, f);
  return ss.str();
}

std::string PLLUnrootedTree::getNewickString(UnodePrinter f,
      corax_unode_t *root, 
      bool rooted)
{
  if (!root) {
    root = getAnyInnerNode();
  }
  if (rooted) {
    return getRootedNewickString(root, f);
  } else {
    std::stringstream ss;
    ss << "(";
    printAux(root->back, ss, f);
    auto temp = root->next;
    while (temp != root) {
      ss << ",";
      printAux(temp->back, ss, f);
      temp = temp->next;
    }
    ss << ");";
    return ss.str();
  }
}
  
std::string PLLUnrootedTree::getRootedNewickString(corax_unode_t *root, 
      UnodePrinter f)
{
  std::stringstream ss;
  auto rootLength = root->length;
  auto rootBackLength = root->back->length;
  assert(rootLength == rootBackLength);
  // temporarily divide the virtual root BL by two, to place
  // the root at the middle of this branch in the rooted newick
  root->length /= 2.0;
  root->back->length /= 2.0;
  ss << "(";
  printAux(root, ss, f);
  ss << ",";
  printAux(root->back, ss, f);
  ss << ");";
  root->length = rootLength;
  root->back->length = rootBackLength;
  return ss.str();
}

std::unordered_set<std::string> PLLUnrootedTree::getLabels() const
{
  std::unordered_set<std::string> res;
  for (auto node:  getLeaves()) {
    if (node->label) {
      res.insert(node->label);
    }
  }
  return res;
}

static corax_unode_t *searchForSet(corax_unode_t *node, 
    std::unordered_set<std::string> &currentNodeSet,
    const std::unordered_set<std::string> &set)
{
  if (!node->next) {
    if (node->label) {
      currentNodeSet.insert(std::string(node->label));
    }
  } else {
    auto left = node->next->back;
    auto right = node->next->next->back;
    std::unordered_set<std::string> rightNodeSet;
    auto res1 = searchForSet(left, currentNodeSet, set);
    if (res1) {
      return res1;
    }
    auto res2 = searchForSet(right, rightNodeSet, set);
    if (res2) {
      return res2;
    }
    for (auto &elem: rightNodeSet) {
      currentNodeSet.insert(elem);
    }
  }
  return (currentNodeSet == set) ? node : nullptr;
}


corax_unode_t *PLLUnrootedTree::getVirtualRoot(PLLRootedTree &referenceTree)
{
  std::unordered_set<std::string> leftRLeaves;
  PLLRootedTree::getLeafLabelsUnder(referenceTree.getRoot()->left, leftRLeaves);
  // find a leaf that is not in leftLeaves
  corax_unode_t *virtualRoot = nullptr;
  for (auto leaf: getLeaves()) {
    if (leftRLeaves.find(std::string(leaf->label)) == leftRLeaves.end()) {
      virtualRoot = leaf->back;
      break;
    }
  }
  std::unordered_set<std::string> temp;
  return searchForSet(virtualRoot, temp, leftRLeaves);
}

static size_t leafHash(corax_unode_t *leaf) {
  assert(leaf);
  std::hash<std::string> hash_fn;
  return hash_fn(std::string(leaf->label));
}

static size_t getTreeHashRec(corax_unode_t *node, size_t i) {
  assert(node);
  if (i == 0) 
    i = 1;
  if (!node->next) {
    return leafHash(node);
  }
  auto hash1 = getTreeHashRec(node->next->back, i + 1);
  auto hash2 = getTreeHashRec(node->next->next->back, i + 1);
  std::hash<size_t> hash_fn;
  auto m = std::min(hash1, hash2);
  auto M = std::max(hash1, hash2);
  return hash_fn(m * i + M);

}

static corax_unode_t *findMinimumHashLeafRec(corax_unode_t * root, size_t &hashValue)
{
  assert(root);
  if (!root->next) {
    hashValue = leafHash(root);
    return root;
  }
  size_t hash1, hash2;
  auto temp = root->next;
  auto min1 = findMinimumHashLeafRec(temp->back, hash1);
  temp = temp->next;
  while (temp != root) {
    auto min2 = findMinimumHashLeafRec(temp->back, hash2);
    if (hash2 < hash1) {
      hash1 = hash2;
      min1 = min2;
    }
    temp = temp->next;
  }
  hashValue = hash1;
  return min1;
}

static corax_unode_t *findMinimumHashLeaf(corax_unode_t * root) 
{
  assert(root);
  auto n1 = root;
  auto n2 = root->back;
  size_t hash1 = 0;
  size_t hash2 = 0;
  auto min1 = findMinimumHashLeafRec(n1, hash1);
  auto min2 = findMinimumHashLeafRec(n2, hash2);
  if (hash1 < hash2) {
    return min1;
  } else {
    return min2;
  }
}

size_t PLLUnrootedTree::getUnrootedTreeHash() const
{
  auto minHashLeaf = findMinimumHashLeaf(getAnyInnerNode());
  auto res = getTreeHashRec(minHashLeaf, 0) + getTreeHashRec(minHashLeaf->back, 0);
  return res;
}

size_t PLLUnrootedTree::getRootedTreeHash(corax_unode_t *root) const
{
  auto res = getTreeHashRec(root, 0) + getTreeHashRec(root->back, 0);
  return res;
}



/**
 *  Recursively fill the vector orderedChildren, such that
 *  orderedChildren[node->node_index] is the list of the children
 *  of node ordered by the smallest leaf label under each child.1
 */
static std::string orderChildren(corax_unode_t *node,
    std::vector<std::vector<corax_unode_t *> > &orderedChildren)
{
  if (!node->next) {
    return node->label;
  }
  assert(orderedChildren.size() > node->node_index);
  auto &myChildren = orderedChildren[node->node_index];
  auto temp = node->next;
  std::vector<std::pair<std::string, corax_unode_t *> > toSort;
  while (temp != node) {
    auto firstLabel = orderChildren(temp->back, orderedChildren);
    toSort.push_back({firstLabel, temp->back});
    temp = temp->next;
  }
  std::sort(toSort.begin(), toSort.end());
  for (const auto &pair: toSort) {
    assert(pair.second);
    myChildren.push_back(pair.second);
  }
  return toSort[0].first;
}

/**
 *  returns true if the nodes node1 and node2 are isomorphic
 *  leftFirst1 and leftFirst2 must have been filled with
 *  the orderChildren function
 */
static bool areIsomorphicAux(corax_unode_t *node1,
    corax_unode_t *node2,
    const std::vector< std::vector< corax_unode_t *> > &orderedChildren1,
    const std::vector< std::vector< corax_unode_t *> > &orderedChildren2)
{
  assert(node1);
  assert(node2);
  if (!node1->next || !node2->next) {
    // at least one is a leaf
    if (node1->next || node2->next) {
      // only one is a leaf 
      return false;
    }
    return strcmp(node1->label, node2->label) == 0;
  }
  // both are internal nodes
  if (orderedChildren1.size() != orderedChildren2.size()) {
    return false;
  }
  assert(orderedChildren1.size());
  for (unsigned int i = 0; i < orderedChildren1[node1->node_index].size(); ++i) {
    if (!areIsomorphicAux(
        orderedChildren1[node1->node_index][i],
        orderedChildren2[node2->node_index][i],
        orderedChildren1,
        orderedChildren2)) {
      return false;
    }
  }
  return true;
}

bool PLLUnrootedTree::areIsomorphic(const PLLUnrootedTree &t1,
    const PLLUnrootedTree &t2)
{
  if (t1.getNodeNumber() != t2.getNodeNumber()) {
    return false;
  }
  auto N = t1.getDirectedNodeNumber();
  auto startingLeaf1 = findMinimumHashLeaf(t1.getAnyInnerNode())->back;
  auto startingLeaf2 = findMinimumHashLeaf(t2.getAnyInnerNode())->back;
  std::vector< std::vector< corax_unode_t *> > orderedChildren1(N);
  std::vector< std::vector< corax_unode_t *> > orderedChildren2(N);
  orderChildren(startingLeaf1, orderedChildren1);
  orderChildren(startingLeaf2, orderedChildren2);
  return areIsomorphicAux(startingLeaf1, 
      startingLeaf2, 
      orderedChildren1, 
      orderedChildren2);
}

bool PLLUnrootedTree::isBinary() const
{
  for (auto node: getInnerNodes()) {
    assert(node->next);
    if (node->next->next->next != node) {
      return false;
    }
  }
  return true;
}
  
corax_unode_t *PLLUnrootedTree::findLeaf(const std::string &label)
{
  for (auto leaf: getLeaves()) {
    if (label == leaf->label) {
      return leaf;
    }
  }
  return nullptr;
}

void PLLUnrootedTree::ensureUniqueLabels() 
{
  auto labels = getLabels();
  unsigned int i = 0;
  std::string prefix("s");
  for (auto node: getPostOrderNodes()) {
    if (node->next) {
      std::string newLabel;
      if (node->label) {
        newLabel = std::string(node->label);
      }
      while (labels.find(newLabel) != labels.end() || newLabel.size() == 0) {
        newLabel = prefix + std::to_string(i++);
      }
      free(node->label);
      node->label = static_cast<char*>(malloc(sizeof(char) * (newLabel.size() + 1)));
      std::strcpy(node->label, newLabel.c_str());
      labels.insert(newLabel);
    }
  }
}
  
bool PLLUnrootedTree::hasUniqueLeafLabels() const
{
  std::unordered_set<std::string> labels;
  for (auto node: getLeaves()) {
    std::string label(node->label);
    if (labels.end() != labels.find(label)) {
      return false;
    }
    labels.insert(label);
  }
  return true;
}
  



