#ifndef TRACKRECO_TPCSEEDTRACKMAP_H
#define TRACKRECO_TPCSEEDTRACKMAP_H

#include <trackbase/TrkrDefs.h>


#include <map>
#include <vector>
#include <memory>

class TpcSeedTrackMap
{
 public:

  TpcSeedTrackMap(){};

  std::multimap<unsigned int, unsigned int> SeedTrackMap;

};

#endif
