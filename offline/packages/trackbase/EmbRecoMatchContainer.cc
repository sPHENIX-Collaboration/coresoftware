#include "EmbRecoMatchContainer.h"

namespace {
  EmbRecoMatchContainer::Vector dummy_vec;
  std::vector<unsigned short> dummy_short_vec;

  EmbRecoMatchContainer::Map_nMatches dummy_map;
}

EmbRecoMatchContainer::ConstRange EmbRecoMatchContainer::getMatchedRange() const
{ return make_pair(dummy_vec.cbegin(), dummy_vec.cend()); }

EmbRecoMatchContainer::Vector&    EmbRecoMatchContainer::getMatches() 
{ return dummy_vec; }

std::vector<unsigned short>& EmbRecoMatchContainer::ids_TruthUnmatched() // id's of the TrkrTruthTrack's that are not matched
{ return dummy_short_vec; }

std::map<unsigned short, unsigned short>& EmbRecoMatchContainer::map_nRecoPerTruth() 
{ return dummy_map; }

std::map<unsigned short, unsigned short>& EmbRecoMatchContainer::map_nTruthPerReco() 
{ return dummy_map; }
