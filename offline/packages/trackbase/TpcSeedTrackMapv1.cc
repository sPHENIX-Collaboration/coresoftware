#include "TpcSeedTrackMapv1.h"

#include <TSystem.h>
#include <iostream>

void TpcSeedTrackMapv1::Reset()
{
  SeedTrackMap.clear();
}

void TpcSeedTrackMapv1::addAssoc(unsigned int tpc_key, unsigned int track_key)
{
  SeedTrackMap.insert(std::make_pair(tpc_key, track_key));
}

  // get map entries for one TPC seed
TpcSeedTrackMapv1::ConstRange TpcSeedTrackMapv1::getAssocTracks(unsigned int tpc_key)
  {
    return std::make_pair( SeedTrackMap.lower_bound(tpc_key), SeedTrackMap.upper_bound(tpc_key) ); 
  }

  // get map entries for all TPC seed
TpcSeedTrackMapv1::ConstRange TpcSeedTrackMapv1::getAll()
  {
    return std::make_pair( SeedTrackMap.begin(), SeedTrackMap.end() ); 
  }

unsigned int TpcSeedTrackMapv1::size()
  {
    return SeedTrackMap.size();
  }
