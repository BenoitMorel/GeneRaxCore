#include "UFBoot.hpp"
#include <limits>
#include <maths/Random.hpp>
#include <parallelization/ParallelContext.hpp>
#include <cassert>
#include <IO/Logger.hpp>

Bootstrap::Bootstrap(unsigned int samples)
{
  assert(ParallelContext::isRandConsistent());
  // we generate a subsampling of the elements
  // samples is the number of samples local to the current core
  // we need to subsample over the total number of samples over
  // all cores
  auto totalSamples = samples;
  ParallelContext::sumUInt(totalSamples);
  std::vector<unsigned int> perCoreSamples;
  ParallelContext::allGatherUInt(samples, perCoreSamples);
  unsigned int begin  = 0;
  for (unsigned int i = 0; i < ParallelContext::getRank(); ++i) {
    begin += perCoreSamples[i];  
  }
  auto end = begin + samples;
  for (unsigned int i = 0; i < totalSamples; ++i) {
    unsigned int v = Random::getInt(0, totalSamples - 1);
    if (v >= begin && v < end) {
      indices.push_back(v - begin);
    }
  }
  unsigned int totalSize = indices.size();
  ParallelContext::sumUInt(totalSize);
  assert(totalSize == totalSamples);
}

double Bootstrap::evaluate(const std::vector<double> &likelihoods) const
{
  auto ll = 0.0;
  for (auto i: indices) {
    ll += likelihoods[i];
  }
  ParallelContext::sumDouble(ll);
  return ll;
}


RootBoot::RootBoot(unsigned int samples):
  bootstrap(samples)
{
  reset();
}

void RootBoot::testRoot(const std::vector<double> &values, unsigned int id)
{
  double ll = bootstrap.evaluate(values);
  if (ll > bestLL) {
    bestId = id;
    bestLL = ll;
  }
}

void RootBoot::reset()
{
  bestId = 0;
  bestLL = std::numeric_limits<double>::lowest();
}


PerBranchBoot::PerBranchBoot(unsigned int elements, unsigned int branches):
  _bootstrap(elements),
  _bestLLs(branches),
  _ok(branches)
{
  reset();
}

void PerBranchBoot::test(const std::vector<double> &values, 
    const std::vector<unsigned int> &branches,
    bool isReferenceTree)
{
  double ll = _bootstrap.evaluate(values);
  for (auto branch: branches) {
    
    if (ll > _bestLLs[branch]) {
      //Logger::info << "ko " << ll << " " <<  _bestLLs[branch] << std::endl;
      _bestLLs[branch] = ll;
      _ok[branch] = isReferenceTree;
    } else {
      //Logger::info << "ok " << ll << " " <<  _bestLLs[branch] << std::endl;

    }
  }
}

void PerBranchBoot::reset() 
{
  std::fill(_bestLLs.begin(), _bestLLs.end(), std::numeric_limits<double>::lowest());
  std::fill(_ok.begin(), _ok.end(), true);
}

PerBranchKH::PerBranchKH(unsigned int elements, unsigned int branches, unsigned int bootstraps):
  _perBootstrapRefLL(bootstraps, std::numeric_limits<double>::lowest()),
  _oks(branches, bootstraps)
{
  for (unsigned int i = 0; i < bootstraps; ++i) {
    _bootstraps.push_back(Bootstrap(elements));
  }
  
}

void PerBranchKH::test(const std::vector<double> &values, 
      const std::vector<unsigned int> &branches)
{
  double ll2 = std::accumulate(values.begin(), values.end(), 0.0);
  ParallelContext::sumDouble(ll2);
  if (_refLL - ll2 < -1e-3) {
    Logger::info << "ERROR _refLL - ll2 < -1e-3" << std::endl;
    Logger::info << _refLL << " " << ll2 << std::endl;
  }
  if (_refLL < ll2) {
    ll2 = _refLL;
  }
  double averageDelta = 0.0;
  std::vector<double> deltas(_bootstraps.size(), 0.0);
  for (unsigned int i = 0; i < _bootstraps.size(); ++i) {
    deltas[i] = _perBootstrapRefLL[i] - _bootstraps[i].evaluate(values);
    averageDelta += deltas[i];
  }
  averageDelta /= static_cast<double>(_bootstraps.size());
  unsigned int oks = _bootstraps.size();;
  for (unsigned int i = 0; i < _bootstraps.size(); ++i) {
    if (deltas[i] - averageDelta > _refLL - ll2) {
      oks--;
    }
  }
  for (auto b: branches) {
    _oks[b] = std::min(_oks[b], oks);
  }
}
void PerBranchKH::newML(const std::vector<double> &values)
{
  _refLL = std::accumulate(values.begin(), values.end(), 0.0);
  ParallelContext::sumDouble(_refLL);
  for (unsigned int i = 0; i < _bootstraps.size(); ++i) {
    _perBootstrapRefLL[i] = _bootstraps[i].evaluate(values);
  }
}

void PerBranchKH::newMLTree(const std::vector<double> &values)
{
  newML(values);
  std::fill(_oks.begin(), _oks.end(), _bootstraps.size());
}


