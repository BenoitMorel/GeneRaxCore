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

int main()
{
  // the difference between t1 and t2 is the resolution of (a,b,c)
  std::string t1("(((a,b),c), (d,e), (f,g));");
  std::string t1same("(((b,a),c), (d,e), (g,f));");
  std::string t2("(((a,c),b), (d,e), (f,g));");
  // this should be the strict consensus
  std::string strict("((a,b,c), (d,e), (f,g));");

  PLLUnrootedTree tree1(t1, false);
  PLLUnrootedTree tree2(t2, false);
  assert(!(t1 == t2));
  PLLUnrootedTree treeStrict(strict, false);
  

  {
    auto trees1 = getList(t1, t1same, t2, 20, 20, 20);
    assert(trees1.size() == 60);
    PLLUnrootedTree consensus50(
        PLLUnrootedTree::buildConsensusTree(trees1, 0.5), false);
    assert(consensus50 == tree1); 
    PLLUnrootedTree consensus90(
        PLLUnrootedTree::buildConsensusTree(trees1, 0.9), false);
    Logger::info << consensus90.getNewickString() << std::endl;
    assert(consensus90 == treeStrict);
  }



  return 0;
}
