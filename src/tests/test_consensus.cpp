#include <iostream>
#include <trees/PLLUnrootedTree.hpp>
#include <trees/PLLRootedTree.hpp>

std::vector<std::string> getList(const std::string &t1,
    const std::string &t2,
    const std::string &t3,
    unsigned int n1, 
    unsigned int n2,
    unsigned int n3)
{
  std::vector<std::string> res;
  for (unsigned int i = 0; i < n1; ++i) {
    res.push_back(t1);
  }
  for (unsigned int i = 0; i < n2; ++i) {
    res.push_back(t2);
  }
  for (unsigned int i = 0; i < n3; ++i) {
    res.push_back(t3);
  }
  return res;
}

void testUnrootedConsensus() {
  // the difference between t1 and t2 is the resolution of (a,b,c)
  std::string t1("(((a,b),c), (d,e), (f,g));");
  std::string t1same("(((b,a),c), (d,e), (g,f));");
  std::string t2("(((a,c),b), (d,e), (f,g));");
  // this should be the strict consensus
  std::string strict("((a,b,c), (d,e), (f,g));");

  PLLUnrootedTree tree1(t1, false);
  PLLUnrootedTree tree2(t2, false);
  assert(!(tree1 == tree2));
  PLLUnrootedTree treeStrict(strict, false);
  
  auto trees1 = getList(t1, t1same, t2, 20, 20, 20);
  assert(trees1.size() == 60);
  PLLUnrootedTree consensus50(
      PLLUnrootedTree::buildConsensusTree(trees1, 0.5), false);
  assert(consensus50 == tree1); 
  PLLUnrootedTree consensus90(
      PLLUnrootedTree::buildConsensusTree(trees1, 0.9), false);
  assert(consensus90 == treeStrict);
}

void testRootedConsensus() 
{
  std::string t1 = "((d,(c,(b,a))), (e, (f,g)));";
  std::string t2 = "((e,(d,(b,(a,c)))),(f,g));";
  PLLRootedTree tree1(t1, false);
  PLLRootedTree tree2(t2, false);
  assert(!(tree1 == tree2));

  std::vector<std::string> newicks;
  for (unsigned int i = 0; i < 66; ++i) {
    newicks.push_back(t1);
  }
  for (unsigned int i = 0; i < 34; ++i) {
    newicks.push_back(t2);
  }
  auto consensusStr50 = PLLRootedTree::buildConsensusTree(newicks, 0.5);
  PLLUnrootedTree realConsensus50("((d,(c,(b,a))), (e, (f,g)), o);", false);
  auto consensus50 = PLLUnrootedTree::buildWithOutgroupFromRootedNewick(consensusStr50, "o");
  assert(*consensus50 == realConsensus50);
  assert(!(*consensus50 == tree2));
  
  auto consensusStr90 = PLLRootedTree::buildConsensusTree(newicks, 0.9);
  PLLUnrootedTree realConsensus90("((d,(c,b,a)), e, (f,g), o);", false);
  auto consensus90 = PLLUnrootedTree::buildWithOutgroupFromRootedNewick(consensusStr90, "o");
  assert(*consensus90 == realConsensus90);
  assert(!(*consensus50 == *consensus90));
  


}

int main()
{
  testUnrootedConsensus();
  testRootedConsensus();
  return 0;
}
