#include <iostream>
#include <trees/PLLUnrootedTree.hpp>
#include <trees/PLLRootedTree.hpp>

int main()
{
  std::string rootedStr("((left1,(left2, (left3, left4))),((right1,right2),(right3,right4)));");
  std::string unrootedStr("((left1,(left2, (left3, left4))),((right1,right2),(right3,right4)), outgroup);");
  PLLRootedTree rooted1(rootedStr, false);
  PLLUnrootedTree unrooted1(unrootedStr, false);
  auto unrooted2 = PLLUnrootedTree::buildWithOutgroupFromRooted(rooted1, "outgroup");
  assert(unrooted2 && (*unrooted2 == unrooted1));
 
  auto outgroup = unrooted2->findLeaf("outgroup");
  auto rooted2 = PLLRootedTree::buildFromOutgroup(outgroup);
  assert(rooted2 && *rooted2 == rooted1);
  return 0;
}
