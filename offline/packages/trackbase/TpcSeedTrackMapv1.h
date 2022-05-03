#ifndef TRACKRECO_TPCSEEDTRACKMAPV1_H
#define TRACKRECO_TPCSEEDTRACKMAPV1_H

#include  "TpcSeedTrackMap.h"
#include "TrkrDefs.h"

#include <map>
#include <vector>
#include <memory>

class TpcSeedTrackMapv1 : public TpcSeedTrackMap
{
public:

   using Map = std::multimap<unsigned int, unsigned int>;
   using ConstIterator = Map::const_iterator;
   using ConstRange = std::pair<Map::const_iterator, Map::const_iterator>;

  //! ctor
  TpcSeedTrackMapv1() = default;

  void Reset() override;

  void addAssoc(unsigned int tpc_key, unsigned int track_key) override;

  // get map entries for one TPC seed
  ConstRange getAssocTracks(unsigned int) override;

  ConstRange getAll() override;

  unsigned int size() override;

private:

  Map SeedTrackMap;


 ClassDefOverride(TpcSeedTrackMapv1, 1);

};

#endif
