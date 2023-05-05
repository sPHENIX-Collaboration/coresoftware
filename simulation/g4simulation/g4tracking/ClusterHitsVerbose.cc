#include "ClusterHitsVerbose.h"

namespace {
  ClusterHitsVerbose::Vector dummy_vec;
  std::vector<unsigned short> dummy_short_vec;

  ClusterHitsVerbose::Map_nMatches dummy_map;
}

ClusterHitsVerbose::ConstRange ClusterHitsVerbose::getMatchedRange() const
{ return make_pair(dummy_vec.cbegin(), dummy_vec.cend()); }

ClusterHitsVerbose::Vector&    ClusterHitsVerbose::getMatches() 
{ return dummy_vec; }

std::vector<unsigned short>& ClusterHitsVerbose::ids_TruthUnmatched() // id's of the TrkrTruthTrack's that are not matched
{ return dummy_short_vec; }

std::map<unsigned short, unsigned short>& ClusterHitsVerbose::map_nRecoPerTruth() 
{ return dummy_map; }

std::map<unsigned short, unsigned short>& ClusterHitsVerbose::map_nTruthPerReco() 
{ return dummy_map; }
