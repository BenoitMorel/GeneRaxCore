#include "util/Scenario.hpp"
#include <IO/Logger.hpp>
#include <IO/FileSystem.hpp>
#include <IO/ReconciliationWriter.hpp>
#include <IO/ParallelOfstream.hpp>
#include <vector>
#include <string>
#include <map>
#include <sstream>

const char *Scenario::eventNames[]  = {"S", "SL", "D", "T", "TL", "L", "Leaf", "Invalid"};

const unsigned int Scenario::EVENT_TYPE_NUMBER = 6;

struct TransferPair {
  double count;
  unsigned int id1;
  unsigned int id2;
  TransferPair(double count, unsigned int id1, unsigned int id2):
    count(count), id1(id1), id2(id2) {}
  bool operator < (const TransferPair& p) const
  {
    if (count != p.count) {
      return count < p.count;
    } else if (id1 != p.id1) {
      return id1 < p.id1;
    } else {
      return id2 < p.id2;
    }
  }
};
  
void PerSpeciesEvents::parallelSum() 
{
  for (auto &speciesEvents: events) 
  {
    ParallelContext::sumUInt(speciesEvents.LeafCount);
    ParallelContext::sumUInt(speciesEvents.DCount);
    ParallelContext::sumUInt(speciesEvents.SCount);
    ParallelContext::sumUInt(speciesEvents.SLCount);
    ParallelContext::sumUInt(speciesEvents.TCount);
    ParallelContext::sumUInt(speciesEvents.TLCount);
  }
}

void Scenario::addEvent(ReconciliationEventType type, 
    unsigned int geneNode, 
    unsigned int speciesNode, 
    unsigned int destSpeciesNode) 
{
  
  addTransfer(type, geneNode, speciesNode, destSpeciesNode);
}
  
void Scenario::addTransfer(ReconciliationEventType type, 
  unsigned int geneNode, 
  unsigned int speciesNode, 
  unsigned int destSpeciesNode)
{
  Event event;
  event.type = type;
  event.geneNode = geneNode;
  event.speciesNode = speciesNode;
  event.destSpeciesNode = destSpeciesNode;
  addEvent(event);
}
  
void Scenario::addEvent(const Event &event)
{
  _events.push_back(event);
  assert(static_cast<int>(event.type) >= 0);
  _eventsCount[static_cast<unsigned int>(event.type)] ++;
  if (_geneIdToEvents.size() <= static_cast<size_t>(event.geneNode)) {
    _geneIdToEvents.resize(event.geneNode + 1);
  }
  _geneIdToEvents[event.geneNode].push_back(event);

}

void Scenario::saveEventsCounts(const std::string &filename, bool masterRankOnly) {
  ParallelOfstream os(filename, masterRankOnly);
  for (unsigned int i = 0; i < static_cast<unsigned int>(ReconciliationEventType::EVENT_Invalid); ++i) {
    os << eventNames[i] << ":" << _eventsCount[i] << std::endl;
  }
}


std::vector<unsigned int> Scenario::getPerSpeciesCopies() const
{
  std::vector<unsigned int> copies(_rootedTree->getNodeNumber(), 0);
  for (const auto &event: _events) {
    if (event.type == ReconciliationEventType::EVENT_S
        || event.type == ReconciliationEventType::EVENT_SL
        || event.type == ReconciliationEventType::EVENT_None) {
      copies[event.speciesNode]++;
    }
  }
  return copies;
}
  

void Scenario::gatherReconciliationStatistics(PerSpeciesEvents &perSpeciesEvents) const
{
  for (auto &event: _events) {
    auto e = event.speciesNode;
    auto &speciesEvents = perSpeciesEvents.events[e];
    switch (event.type) {
      case ReconciliationEventType::EVENT_S: 
        speciesEvents.SCount++;
        break;
      case ReconciliationEventType::EVENT_None:
        speciesEvents.LeafCount++;
        break;
      case ReconciliationEventType::EVENT_SL: 
        speciesEvents.SLCount++;
        break;
      case ReconciliationEventType::EVENT_D:
        speciesEvents.DCount++;
        break;
      case ReconciliationEventType::EVENT_T:
        speciesEvents.TCount++;
        break;
      case ReconciliationEventType::EVENT_TL:
        speciesEvents.TLCount++;
        break;
      default:
        assert(false);
    }
  }
}


static void dumpSpeciesEventCountVector(ParallelOfstream &os,
    const std::vector<double> &eventCounts,
    const std::vector<std::string> &indexToLabel)
{
  os << "species_label, speciations, duplications, losses, transfers, presence, origination" << std::endl;
  auto N = indexToLabel.size();
  assert(Scenario::EVENT_TYPE_NUMBER * N == eventCounts.size());
  for (unsigned int i = 0; i < N; ++i) {
    bool print = false;
    for (unsigned int e = 0; e < Scenario::EVENT_TYPE_NUMBER; ++e) {
      print |= (eventCounts[i * Scenario::EVENT_TYPE_NUMBER + e] > 0.0);
    }
    if (print) {
      os << indexToLabel[i];
      for (unsigned int e = 0; e < Scenario::EVENT_TYPE_NUMBER; ++e) {
        os << ", " <<  eventCounts[i * Scenario::EVENT_TYPE_NUMBER + e];
      }
      os << "\n";
    }
  }
}

static void dumpSpeciesToEventCount(ParallelOfstream &os,
    const std::unordered_map<std::string, std::vector<double> > &speciesToEventCount)
{
  os << "species_label, speciations, duplications, losses, transfers, presence, origination" << std::endl;
  std::vector<double> defaultCount(6, 0.0);
  for (auto &it: speciesToEventCount) {
    if (defaultCount == it.second) {
      // do not write species without any event
      continue;
    }
    assert(it.second.size() == 6);
    os << it.first;
    for (auto v: it.second) {
      os << ", " << v;
    }
    os << "\n";
  } 
}

void Scenario::savePerSpeciesEventsCounts(const std::string &filename, bool masterRankOnly) {
  
  ParallelOfstream os(filename, masterRankOnly);
  std::unordered_map<std::string, std::vector<double> > speciesToEventCount;
  std::vector<double> defaultCount(6, 0);
  for (unsigned int e = 0; e < _speciesTree->tip_count + _speciesTree->inner_count; ++e) {
    assert(_speciesTree->nodes[e]->label);
    speciesToEventCount.insert({std::string(_speciesTree->nodes[e]->label), defaultCount});
  }
  // 0: speciation (or numbero of gene copies)
  // 1: duplication
  // 2: loss:
  // 3: transfer
  // 4: presence (1 or 0)
  // 5: origination
  for (auto &event: _events) {
    auto &eventCount = speciesToEventCount[_speciesTree->nodes[event.speciesNode]->label];
    switch (event.type) {
      case ReconciliationEventType::EVENT_S: 
      case ReconciliationEventType::EVENT_None:
        eventCount[0]++;
        eventCount[4] = 1;
        break;
      case ReconciliationEventType::EVENT_SL: 
        eventCount[0]++;
        eventCount[4] = 1;
        // count the loss
        speciesToEventCount[event.pllLostSpeciesNode->label][2]++;
        break;
      case ReconciliationEventType::EVENT_D:
        eventCount[1]++;
        break;
      case ReconciliationEventType::EVENT_T:
        eventCount[3]++;
        break;
      case ReconciliationEventType::EVENT_TL:
        eventCount[2]++;
        eventCount[3]++;
        break;
      case ReconciliationEventType::EVENT_L: case ReconciliationEventType::EVENT_Invalid:
        break;
    }
  }
  auto origin = getOriginationSpecies();
  assert(origin && origin->label);
  speciesToEventCount[origin->label][5] = 1.0;
  dumpSpeciesToEventCount(os, speciesToEventCount);
}
  
corax_rnode_t *Scenario::getOriginationSpecies() const
{
  auto &rootEvents = _geneIdToEvents[_virtualRootIndex];
  assert(rootEvents.size());
  return _speciesTree->nodes[rootEvents[0].speciesNode];
}

void Scenario::mergeTransfers(const PLLRootedTree &speciesTree,
    const std::string &filename,
    const std::vector<std::string> &filenames,
    bool parallel,
    bool normalize)
{
  ParallelOfstream os(filename, parallel);
  const auto labelToId = speciesTree.getDeterministicLabelToId();
  const auto idToLabel = speciesTree.getDeterministicIdToLabel();
  const unsigned int N = labelToId.size();
  const VectorDouble zeros(N, 0.0);
  auto countMatrix = MatrixDouble(N, zeros);
  for (const auto &f: filenames) {
    std::ifstream is(f);
    std::string line;
    while (std::getline(is, line)) {
      if (line[0] == '#') {
        continue;
      }
      std::istringstream iss(line);
      std::string sp1;
      std::string sp2;
      double count = 0.0;
      iss >> sp1 >> sp2 >> count;
      countMatrix[labelToId.at(sp1)][labelToId.at(sp2)] += count;
    }
  }
  unsigned int subfileCount = filenames.size();
  if (parallel) {
    ParallelContext::sumUInt(subfileCount);
  }
  std::vector<TransferPair> transfers;
  for (unsigned int i = 0; i < N; ++i) {
    for (unsigned int j = 0; j < N; ++j) {
      if (parallel) {
        ParallelContext::sumDouble(countMatrix[i][j]);
      }
      if (0.0 == countMatrix[i][j]) {
        continue;
      }
      double transferCount = countMatrix[i][j];
      if (normalize) {
        transferCount /= static_cast<double>(subfileCount);
      }
      transfers.push_back(TransferPair(transferCount, i, j));
    }
  }
  std::sort(transfers.rbegin(), transfers.rend());
  for (const auto &t: transfers) {
    os << idToLabel[t.id1] << " " << idToLabel[t.id2] << " " << t.count << std::endl;

  }
}


void Scenario::mergePerSpeciesEventCounts(const PLLRootedTree &speciesTree,
    const std::string &filename,
    const std::vector<std::string> &filenames,
    bool parallel,
    bool normalize)
{
  ParallelOfstream os(filename, parallel);
  auto speciesLabelToIndex = speciesTree.getDeterministicLabelToId();
  auto speciesIndexToLabel = speciesTree.getDeterministicIdToLabel();
  auto N = speciesLabelToIndex.size();
  std::vector<double> eventCounts(EVENT_TYPE_NUMBER * N, 0.0); 
  for (const auto &subfile: filenames) {
    std::ifstream is(subfile);
    std::string line;
    std::getline(is, line);
    while (std::getline(is, line)) {
      if (line[0] == '#') {
        continue;
      }
      std::istringstream iss(line);
      std::string species;
      iss >> species;
      species.pop_back(); // remove comma
      auto speciesIndex = speciesLabelToIndex[species];
      auto begin = speciesIndex * EVENT_TYPE_NUMBER;
      for (unsigned int i = begin; i < begin + EVENT_TYPE_NUMBER; ++i) {
        double temp;
        iss >> temp;
        std::string comma;
        iss >> comma;
        eventCounts[i] += temp;
      }
    }
  }
  unsigned int subfileCount = filenames.size(); 
  if (parallel) {
    if (normalize) {
      ParallelContext::sumUInt(subfileCount);
    }
    ParallelContext::sumVectorDouble(eventCounts);
  }
  if (normalize) {
    for (auto &count: eventCounts) {
      count /= static_cast<double>(subfileCount);
    }
  }
  dumpSpeciesEventCountVector(os, eventCounts, speciesIndexToLabel);
}

void Scenario::saveReconciliation(const std::string &filename, ReconciliationFormat format, bool masterRankOnly)
{
  ParallelOfstream os(filename, masterRankOnly);
  saveReconciliation(os, format);
}

void Scenario::saveReconciliation(ParallelOfstream &os, ReconciliationFormat format)
{
  switch (format) {
  case ReconciliationFormat::NHX:
    ReconciliationWriter::saveReconciliationNHX(_speciesTree, 
        _geneRoot, 
        _virtualRootIndex, 
        _geneIdToEvents, 
        os);
    break;
  case ReconciliationFormat::ALE:
    ReconciliationWriter::saveReconciliationALE(_speciesTree, 
        _geneRoot, 
        _virtualRootIndex, 
        _geneIdToEvents, 
        os);
    break;
  case ReconciliationFormat::RecPhyloXML:
    ReconciliationWriter::saveReconciliationRecPhyloXML(_speciesTree, 
        _virtualRootIndex, 
        _geneIdToEvents, 
        os);
    break;
  case ReconciliationFormat::NewickEvents:
    ReconciliationWriter::saveReconciliationNewickEvents( 
        _geneRoot, 
        _virtualRootIndex, 
        _geneIdToEvents, 
        os);
    break;
  }
}
  
void Scenario::countTransfers(const StringToUint &labelToId,
      MatrixUint &count)
{
  char *labelFrom = nullptr;
  char *labelTo = nullptr;
  unsigned int from = 0;
  unsigned int to = 0;
  for (auto &event: _events) {
    switch(event.type) {
    case ReconciliationEventType::EVENT_T:
    case ReconciliationEventType::EVENT_TL:
      labelFrom = _speciesTree->nodes[event.speciesNode]->label;
      labelTo = _speciesTree->nodes[event.destSpeciesNode]->label;
      from = labelToId.at(std::string(labelFrom));
      to = labelToId.at(std::string(labelTo));
      count[from][to]++;
      break;
    case ReconciliationEventType::EVENT_S:
    case ReconciliationEventType::EVENT_SL:
    case ReconciliationEventType::EVENT_L:
    case ReconciliationEventType::EVENT_D:
    case ReconciliationEventType::EVENT_None:
    case ReconciliationEventType::EVENT_Invalid:
      break;
    }
  }
}

void Scenario::countOrigins(const StringToUint &labelToId,
    std::vector<unsigned int> &fromS,
    std::vector<unsigned int> &fromSButL,
    MatrixUint &count)
{
  // origins come from SL, S, T, and TL events
  countTransfers(labelToId, count);
  char *destLabel1 = nullptr;
  char *destLabel2 = nullptr;
  for (auto &event: _events) {
    switch(event.type) {
    case ReconciliationEventType::EVENT_S:
      destLabel1 = _speciesTree->nodes[event.speciesNode]->left->label;
      destLabel2 = _speciesTree->nodes[event.speciesNode]->right->label;
      fromS[labelToId.at(destLabel1)]++;
      fromS[labelToId.at(destLabel2)]++;
      break;
    case ReconciliationEventType::EVENT_SL:
      destLabel1 = event.pllDestSpeciesNode->label; 
      destLabel2 = event.pllLostSpeciesNode->label; 
      assert(destLabel1 != destLabel2);
      fromS[labelToId.at(destLabel1)]++;
      fromSButL[labelToId.at(destLabel2)]++;
      break;
    case ReconciliationEventType::EVENT_T:
    case ReconciliationEventType::EVENT_TL:
    case ReconciliationEventType::EVENT_L:
    case ReconciliationEventType::EVENT_D:
    case ReconciliationEventType::EVENT_None:
    case ReconciliationEventType::EVENT_Invalid:
      break;
    }
  }
}
  
void Scenario::saveTransfers(const std::string &filename, bool masterRankOnly)
{
  ParallelOfstream os(filename, masterRankOnly);
  for (auto &event: _events) {
    if (event.type == ReconciliationEventType::EVENT_T || event.type == ReconciliationEventType::EVENT_TL) {
      os << _speciesTree->nodes[event.speciesNode]->label << " " 
        << _speciesTree->nodes[event.destSpeciesNode]->label << " " << 1  << std::endl;
    }
  }
}
 

double Scenario::countTransfer(const std::string &from, const std::string &to) 
{
  double res = 0.0;
  for (auto &event: _events) {
    if (event.type == ReconciliationEventType::EVENT_T || event.type == ReconciliationEventType::EVENT_TL) {
      if (_speciesTree->nodes[event.speciesNode]->label == from && 
          _speciesTree->nodes[event.destSpeciesNode]->label == to) {
        res += 1.0;
      }
    }
  }
  return res;
}

  
void Scenario::saveLargestOrthoGroup(std::string &filename, bool masterRankOnly) const
{
  ParallelOfstream os(filename, masterRankOnly);
  corax_unode_t virtualRoot;
  virtualRoot.next = _geneRoot;
  virtualRoot.node_index = _virtualRootIndex;
  virtualRoot.label = nullptr;
  virtualRoot.length = 0.0;
  std::unique_ptr<OrthoGroup> orthoGroup(getLargestOrthoGroupRec(&virtualRoot, true));
  for (auto value: *orthoGroup) {
    os << value << std::endl;
  }
  os << "-" << std::endl;

}

void Scenario::saveAllOrthoGroups(std::string &filename, bool masterRankOnly) const
{
  ParallelOfstream os(filename, masterRankOnly);
  corax_unode_t virtualRoot;
  virtualRoot.next = _geneRoot;
  virtualRoot.node_index = _virtualRootIndex;
  virtualRoot.label = nullptr;
  virtualRoot.length = 0.0;
  OrthoGroupPtr currentOrthoGroup = std::make_shared<OrthoGroup>();
  OrthoGroups orthoGroups;
  getAllOrthoGroupRec(&virtualRoot, orthoGroups, currentOrthoGroup, true);
  orthoGroups.push_back(currentOrthoGroup);
  for (auto orthoGroup: orthoGroups) {
    for (auto value: *orthoGroup) {
      os << value << std::endl;
    }
    os << "-" << std::endl;
  }

}



void Scenario::saveTransferPairCountGlobal(PLLRootedTree &speciesTree,
    std::vector< std::shared_ptr<Scenario> > &scenarios,
    const std::string &filename)
{
  const auto labelToId = speciesTree.getDeterministicLabelToId();
  const auto idToLabel = speciesTree.getDeterministicIdToLabel();
  const unsigned int N = labelToId.size();
  const VectorUint zeros(N, 0);
  auto countMatrix = MatrixUint(N, zeros);
  for (auto &scenario: scenarios) {
    scenario->countTransfers(labelToId, countMatrix);
  }
  std::vector<TransferPair> transferPairs;
  for (unsigned int i = 0; i < N; ++i) {
    for (unsigned int j = 0; j < N; ++j) {
      ParallelContext::sumUInt(countMatrix[i][j]);
      if (countMatrix[i][j]) { 
        TransferPair p(countMatrix[i][j], i, j);
        transferPairs.push_back(p);
      }
    }
  }
  unsigned int parentTransfers = 0;
  for (auto speciesNode: speciesTree.getNodes()) {
    std::string fromLabel = speciesNode->label;
    auto from = labelToId.at(fromLabel);
    auto parent = speciesNode;
    while (parent) {
      std::string toLabel = parent->label;
      auto to = labelToId.at(toLabel);
      parentTransfers += countMatrix[from][to];
      parent = parent->parent;
    }
  }
  std::sort(transferPairs.rbegin(), transferPairs.rend());
  ParallelOfstream os(filename);
  for (auto p: transferPairs) {
    os << idToLabel[p.id1] 
      << " " << idToLabel[p.id2] 
      << " " << p.count 
      << std::endl;
  }
}


void Scenario::saveOriginsGlobal(PLLRootedTree &speciesTree,
    std::vector< std::shared_ptr<Scenario> > &scenarios,
    unsigned int samples,
    const std::string &outputDir)
{
  const auto labelToId = speciesTree.getDeterministicLabelToId();
  const auto idToLabel = speciesTree.getDeterministicIdToLabel();
  const unsigned int N = labelToId.size();
  const VectorUint zeros(N, 0);
  auto countMatrix = MatrixUint(N, zeros);
  std::vector<unsigned int> fromS(N, 0);
  std::vector<unsigned int> fromSButL(N, 0);
  for (auto &scenario: scenarios) {
    scenario->countOrigins(labelToId, fromS, fromSButL, countMatrix);
  }
  // iterate over all dest species, and compute their origins
  for (unsigned int j = 0; j < N; ++j) {
    // will stor the origins
    std::vector<TransferPair> transferPairs;
    auto label = idToLabel[j];
    auto output = FileSystem::joinPaths(outputDir, label + ".txt");
    ParallelContext::sumUInt(fromS[j]);
    ParallelContext::sumUInt(fromSButL[j]);
    // iterate over all origins
    for (unsigned int i = 0; i < N; ++i) {
      ParallelContext::sumUInt(countMatrix[i][j]);
      if (countMatrix[i][j]) { 
        TransferPair p(countMatrix[i][j], i, j);
        transferPairs.push_back(p);
      }
    }
    std::sort(transferPairs.rbegin(), transferPairs.rend());
    ParallelOfstream os(output);
    os << "from_parent_surviving, " << double(fromS[j]) / double(samples) << std::endl;
    os << "from_parent_exctinct, " << double(fromSButL[j]) / double(samples) << std::endl;
    for (auto p: transferPairs) {
      os << idToLabel[p.id1] 
        << ", " << double(p.count) / double(samples)
        << std::endl;
    }
  }
}


OrthoGroup *Scenario::getLargestOrthoGroupRec(corax_unode_t *geneNode, bool isVirtualRoot) const
{
  auto &events = _geneIdToEvents[geneNode->node_index];
  for (auto &event: events) {
    if (event.type == ReconciliationEventType::EVENT_TL) {
      return new OrthoGroup();
    }
  }
  if (geneNode->next == nullptr) {
    auto res = new OrthoGroup();
    res->insert(std::string(geneNode->label));
    return res;
  } else {
    auto left = geneNode->next->back;
    auto right = geneNode->next->next->back;
    if (isVirtualRoot) {
      left = geneNode->next;
      right = geneNode->next->back;
    }
    auto leftOrthoGroup = getLargestOrthoGroupRec(left, false);
    auto rightOrthoGroup = getLargestOrthoGroupRec(right, false);
    auto &event = events.back();
    switch (event.type) {
    case ReconciliationEventType::EVENT_S:
      // keep both groups
      for (auto &v: *leftOrthoGroup) {
        rightOrthoGroup->insert(v);
      }
      delete leftOrthoGroup;
      return rightOrthoGroup;
    case ReconciliationEventType::EVENT_T:
      // only keep the non transfered gene
      if (event.rightGeneIndex == left->node_index) {
        delete leftOrthoGroup;
        return  rightOrthoGroup;
      } else if (event.rightGeneIndex == right->node_index) {
        delete rightOrthoGroup;
        return leftOrthoGroup;
      } else {
        assert(false);
      }
    case ReconciliationEventType::EVENT_TL:
    case ReconciliationEventType::EVENT_None:
    case ReconciliationEventType::EVENT_SL:
    case ReconciliationEventType::EVENT_L:
      // we already handled these cases
      // or they should not be the last event
      // attached to this gene node
      assert(false);
      break;
    case ReconciliationEventType::EVENT_Invalid:
      delete leftOrthoGroup;
      delete rightOrthoGroup;
      break;
    case ReconciliationEventType::EVENT_D:
      // keep the bigger group!
      if (leftOrthoGroup->size() < rightOrthoGroup->size()) {
        std::swap(leftOrthoGroup, rightOrthoGroup);
      }
      delete rightOrthoGroup;
      return leftOrthoGroup;
    }
    return new OrthoGroup();
  }
}

static void appendOrtho(OrthoGroupPtr &orthoGroups, 
    const OrthoGroupPtr &orthoGroupsToAppend)
{
  orthoGroups->insert(orthoGroupsToAppend->begin(),
      orthoGroupsToAppend->end());
}

void Scenario::getAllOrthoGroupRec(corax_unode_t *geneNode,
      OrthoGroups &orthoGroups,
      OrthoGroupPtr &currentOrthoGroup,
      bool isVirtualRoot) const
{
  auto &events = _geneIdToEvents[geneNode->node_index];
  bool underTL = false;
  for (auto &event: events) {
    if (event.type == ReconciliationEventType::EVENT_TL) {
      underTL = true;
    }
  }
  if (geneNode->next == nullptr) {
    currentOrthoGroup->insert(std::string(geneNode->label));
  } else {
    auto left = geneNode->next->back;
    auto right = geneNode->next->next->back;
    if (isVirtualRoot) {
      left = geneNode->next;
      right = geneNode->next->back;
    }
    OrthoGroupPtr leftOrthoGroup = std::make_shared<OrthoGroup>();
    OrthoGroupPtr rightOrthoGroup = std::make_shared<OrthoGroup>();
    getAllOrthoGroupRec(left, orthoGroups, leftOrthoGroup, false);
    getAllOrthoGroupRec(right, orthoGroups, rightOrthoGroup, false);
    auto &event = events.back();
    switch (event.type) {
    case ReconciliationEventType::EVENT_S:
      // merge both groups
      appendOrtho(currentOrthoGroup, leftOrthoGroup);
      appendOrtho(currentOrthoGroup, rightOrthoGroup);
      break;
    case ReconciliationEventType::EVENT_T:
      // save the orthoGroup from the transfered gene,
      // and keep filling the orthoGroup on the non transfered gene
      if (event.rightGeneIndex == left->node_index) {
        appendOrtho(currentOrthoGroup, rightOrthoGroup);
        orthoGroups.push_back(leftOrthoGroup);
      } else if (event.rightGeneIndex == right->node_index) {
        appendOrtho(currentOrthoGroup, leftOrthoGroup);
        orthoGroups.push_back(rightOrthoGroup);
      } else {
        assert(false);
      }
      break;
    case ReconciliationEventType::EVENT_TL:
    case ReconciliationEventType::EVENT_None:
    case ReconciliationEventType::EVENT_SL:
    case ReconciliationEventType::EVENT_L:
      // we already handled these cases
      // or they should not be the last event
      // attached to this gene node
      assert(false);
      break;
    case ReconciliationEventType::EVENT_Invalid:
      orthoGroups.push_back(leftOrthoGroup);
      orthoGroups.push_back(rightOrthoGroup);
      break;
    case ReconciliationEventType::EVENT_D:
      // save the biggest group and continue working on the other
      if (leftOrthoGroup->size() < rightOrthoGroup->size()) {
        std::swap(leftOrthoGroup, rightOrthoGroup);
      }
      orthoGroups.push_back(rightOrthoGroup);
      appendOrtho(currentOrthoGroup, leftOrthoGroup); 
      break;
    }
  }
  if (underTL) {
    // we got transfered: separate the current orthogroup
    orthoGroups.push_back(currentOrthoGroup);
    currentOrthoGroup->clear();
  }
}


void Scenario::initBlackList(unsigned int genesNumber, unsigned int speciesNumber)
{
  std::vector<int>  blackListAux(speciesNumber, false);
  _blacklist = std::make_unique<ScenarioBlackList>(genesNumber, blackListAux);
  
}

void Scenario::blackList(unsigned int geneNode, unsigned int speciesNode)
{
  if (_blacklist && _blacklist->size() > geneNode) { // not true for virtual nodes
    (*_blacklist)[geneNode][speciesNode] = true;
  }
}
bool Scenario::isBlacklisted(unsigned int geneNode, unsigned int speciesNode)
{
  if (_blacklist && _blacklist->size() > geneNode) { // not true for virtual nodes
    return (*_blacklist)[geneNode][speciesNode];
  }
  return false;
}

void Scenario::resetBlackList()
{
  for (auto &aux: *_blacklist) {
    for (unsigned int i = 0; i < aux.size(); ++i) {
      aux[i] = false;
    }
  }
}

static corax_unode_t *createNode(unsigned int index)
{
  auto node = new corax_unode_t();
  node->next = nullptr;
  node->back = nullptr;
  node->data = nullptr;
  node->label = nullptr;
  node->node_index = index;
  return node;
} 

static void backLink(corax_unode_t *n1, corax_unode_t *n2)
{
  n1->back = n2;
  n2->back = n1;
}

static void nextLink(corax_unode_t *n1, corax_unode_t *n2, corax_unode_t *n3)
{
  n1->next = n2;
  n2->next = n3;
  n3->next = n1;
}

corax_unode_t *Scenario::generateVirtualGeneRoot()
{
  auto node1 = createNode(0);
  _geneNodeBuffer.push_back(node1);
  return node1;
}


void Scenario::generateGeneChildren(corax_unode_t *geneNode, 
    corax_unode_t *&leftGeneNode,
    corax_unode_t *&rightGeneNode)
{
  if (_geneNodeBuffer.size() == 1) {
    // we are at the root
    leftGeneNode = _geneNodeBuffer.back();
    rightGeneNode = createNode(_geneNodeBuffer.size());
    backLink(leftGeneNode, rightGeneNode);
    _geneNodeBuffer.push_back(rightGeneNode);
    return;
  }
  auto next1 = createNode(geneNode->node_index);
  auto next2 = createNode(geneNode->node_index);
  nextLink(geneNode, next1, next2);
  leftGeneNode = createNode(_geneNodeBuffer.size());
  backLink(next1, leftGeneNode);
  _geneNodeBuffer.push_back(leftGeneNode);
  rightGeneNode = createNode(_geneNodeBuffer.size());
  backLink(next2, rightGeneNode);
  _geneNodeBuffer.push_back(rightGeneNode);
}




