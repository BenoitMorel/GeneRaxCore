#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <set>
#include <maths/bitvector.hpp>
#include <iostream>
#include <util/enums.hpp>

using CID = unsigned int;

using CCPClade = genesis::utils::Bitvector;
using CladeCounts = std::unordered_map<CID, double>;
using SubcladeCounts = std::vector<CladeCounts>;
using OrderedClades = std::set<CCPClade>;

struct CladeSplit {
  CID parent;
  CID left;
  CID right;
  double frequency;
  double deviation;
  double blLeft;
  double blRight;
  CladeSplit():
    frequency(0.0),
    blLeft(0.0),
    blRight(0.0)
  {
  }
};


/**
 *  SubladeBLs[cid][childcid] = vector of branch lengths
 */
// TODO we could replace this by a pair {sum, count} to compute the average!!
using BLs = std::vector<double>;
using CladeBLs = std::unordered_map<unsigned int, BLs>;
using SubcladeBLs = std::vector<CladeBLs>;

using CladeSplits = std::vector<CladeSplit>;
using CladeToCID = std::unordered_map<CCPClade, CID>;
using CIDToClade = std::vector<CCPClade>;
using CIDToLeaf = std::unordered_map<unsigned int, std::string>;
  


/**
 * Conditional clade representation of distribution of gene trees,
 * used to compute conditional clade probabilities.
 * 
 * A clade is represented with a bitvector, in which each bit
 * corresponds to the presence/absence of a given taxon in the calde
 *
 * The ID of a leaf/taxon corresponds to its index in the 
 * bitvectors representing clades
 *
 * The CID is a unique identifier for a clade
 * 
 * A CladeSplit is a split of a clade into two child clades. Its
 * frequency corresponds to how often this clade is split into
 * those two child clades in the distribution of gene trees. 
 *
 * For each non-trivial clade, we store a vector of CladeSplits, 
 * whose frequencies should sum to one.
 *
 */
class ConditionalClades {
 public:
  ConditionalClades() {}
  ConditionalClades(const std::string &inputFile,
      const std::string &likelihoods,
      CCPRooting ccpRooting,
      unsigned int sampleFrequency = 1);
  void printContent() const; 
  void printStats() const;
  unsigned int getCladesNumber() const {return _cladeToCID.size();}
  const CladeSplits &getCladeSplits(CID cid) const {return _allCladeSplits[cid];} 
  bool isLeaf(CID cid) const;
  std::string getLeafLabel(CID cid) const;
  unsigned int getRootsNumber() const;
  unsigned int getInputTreesNumber() const {return _inputTrees;}
  unsigned int getUniqueInputTreesNumber() const {return _uniqueInputTrees;}
  unsigned int getLeafNumber() const {return _CIDToLeaf.size();}
  const std::unordered_map<unsigned int, std::string> &
    getCidToLeaves() {return _CIDToLeaf;}
  bool skip() const {return false;}
  //bool skip() const {return  _uniqueInputTrees == _inputTrees;} 

  void serialize(const std::string &outputFile);
  void unserialize(const std::string &inputFile);
  void buildFromGeneTrees(const std::string &inputFile, 
      const std::string &likelihoods,
      CCPRooting ccpRooting,
      unsigned int sampleFrequency);
  void buildFromALEFormat(const std::string &inputFile, 
        CCPRooting ccpRooting);

  bool madRooting() const {return _ccpRooting == CCPRooting::MAD;}

  void reorderClades(const std::vector<CID> &mappings);

  bool isValid() const {return _isValid;}
private:
  
  unsigned int _inputTrees;
  unsigned int _uniqueInputTrees;
  CCPRooting _ccpRooting;
  std::vector<std::string> _idToLeaf;
  CIDToLeaf _CIDToLeaf;
  CladeToCID _cladeToCID;
  CIDToClade _CIDToClade;
  std::vector<CladeSplits> _allCladeSplits;
  bool _isValid;
private:
  void _fillCCP(SubcladeCounts &subcladeCounts,
      SubcladeBLs &subcladeBLs,
      bool useLikelihoods,
      std::unordered_map<unsigned int, double> *CIDToDeviation = nullptr);
};



