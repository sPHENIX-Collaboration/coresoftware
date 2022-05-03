#ifndef TRACKRECO_TPCSEEDTRACKMAP_H
#define TRACKRECO_TPCSEEDTRACKMAP_H

#include "TrkrDefs.h"

#include  <phool/PHObject.h>

#include <map>
#include <vector>
#include <memory>

class TpcSeedTrackMap : public PHObject
{
public:

   using Map = std::multimap<unsigned int, unsigned int>;
   using ConstIterator = Map::const_iterator;
   using ConstRange = std::pair<Map::const_iterator, Map::const_iterator>;

  void Reset() override;

  virtual void addAssoc(unsigned int tpc_key, unsigned int track_key) = 0;

  virtual ConstRange getAssocTracks(unsigned int) = 0;

  virtual ConstRange getAll() = 0;

  virtual unsigned int size() = 0;

  protected:

  TpcSeedTrackMap() = default;

private:

 ClassDefOverride(TpcSeedTrackMap, 1);

};

#endif
