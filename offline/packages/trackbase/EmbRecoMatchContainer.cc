#include "EmbRecoMatchContainer.h"

namespace {
  EmbRecoMatchContainer::Vector dummy_vec;
  std::vector<unsigned short> dummy_short_vec;
}


EmbRecoMatchContainer::ConstRange EmbRecoMatchContainer::getMatchedRange() const
{ return make_pair(dummy_vec.cbegin(), dummy_vec.cend()); }

EmbRecoMatchContainer::Vector&    EmbRecoMatchContainer::getMatches() 
{ return dummy_vec; }

std::vector<unsigned short>& EmbRecoMatchContainer::ids_Unmatched() // id's of the TrkrTruthTrack's that are not matched
{ return dummy_short_vec; }
