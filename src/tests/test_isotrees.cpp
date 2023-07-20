#include <iostream>
#include <trees/PLLUnrootedTree.hpp>


void testBinaryRef()
{
  PLLUnrootedTree ref("((a,(b,c)), ((d, e), (f, g)), (h, i));", false);
  PLLUnrootedTree ok1("((a,(c,b)), (((d, e), (f, g)), (h, i)));", false);
  PLLUnrootedTree ok2("(((c,b), a), (((d, e), (f, g)), (h, i)));", false);
  PLLUnrootedTree ok3("((a,(b,c)),(i, h), ((d, e), (f, g)));", false);
  // different taxon number
  PLLUnrootedTree ko1("(((d, e), (f, g)), (h, i));", false);
  // different topology
  PLLUnrootedTree ko2("((b,(a,c)), ((d, e), (f, g)), (h, i));", false);
  // this one is multifurcating
  PLLUnrootedTree ko3("((a,b,c), ((d, e), (f, g)), (h, i));", false);
  // different leaf set
  PLLUnrootedTree ko4("((oups,(b,c)), ((d, e), (f, g)), (h, i));", false);

  assert(ref == ok1);
  assert(ref == ok2);
  assert(ref == ok3);
  assert(!(ref == ko1));
  assert(!(ref == ko2));
  assert(!(ref == ko3));
  assert(!(ref == ko4));

}

void testMultifurcatingRef()
{
  PLLUnrootedTree ref("((a,b,c), ((d, e), (f, g), (x,y)), (h, i));", false);
  PLLUnrootedTree ok1("((c,b,a), ((e, d), (g, f), (x,y)), (h, i));", false);
  PLLUnrootedTree ok2("((h,i), (a,b,c), ((d, e), (f, g), (x,y)))plop:0.1;", false);
  // wrong taxon number
  PLLUnrootedTree ko1("((a,b), ((d, e), (f, g), (x,y)), (h, i));", false);
  // different topology
  PLLUnrootedTree ko2("((a,b,c), ((d, e, f), g, (x,y)), (h, i));", false);
  // binary
  PLLUnrootedTree ko3("(((a,b),c), ((d, e), (f, g), (x,y)), (h, i));", false);
  // different leaf set
  PLLUnrootedTree ko4("((coucou,b,c), ((d, e), (f, g), (x,y)), (h, i));", false);
  assert(ref == ok1);
  assert(ref == ok2);
  assert(!(ref == ko1));
  assert(!(ref == ko2));
  assert(!(ref == ko3));
  assert(!(ref == ko4));
}

int main()
{
  testBinaryRef();
  testMultifurcatingRef();
  return 0;
}
