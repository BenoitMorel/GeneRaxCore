#include "GeneRaxCheckpoint.hpp"

#include <IO/ParallelOfstream.hpp>
#include <trees/JointTree.hpp>

GeneRaxCheckpoint::GeneRaxCheckpoint(const std::string checkpointPath):
  path(checkpointPath),
  checkpointExists(false),
  modelParamDone(false)
{
  std::ifstream is(path);
  if (is.fail()) {
    return;
  } 
  std::cout << "checkpoint exists " << std::endl;
  is >> geneTreeNewickStr;
  is >> substModelStr;
  std::string line;
  std::getline(is, line);
  unsigned int dim = 0;
  is >> dim;
  ratesVector = Parameters(dim);
  for (unsigned int i = 0; i < ratesVector.dimensions(); ++i) {
    is >> ratesVector[i];
  } 
  checkpointExists = true;
}
  
void GeneRaxCheckpoint::save(bool masterRankOnly)
{
  assert(geneTreeNewickStr.size());
  assert(substModelStr.size());
  ParallelOfstream os(path, masterRankOnly);
  os << geneTreeNewickStr << std::endl; 
  os << substModelStr << std::endl;
  os << ratesVector.dimensions();
  for (unsigned int i = 0; i < ratesVector.dimensions(); ++i) {
    os << " " << ratesVector[i];
  }
  os << " " << std::endl;
}

