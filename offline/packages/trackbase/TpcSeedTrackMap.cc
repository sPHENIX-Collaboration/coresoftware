#include "TpcSeedTrackMap.h"

#include <TSystem.h>
#include <iostream>

void TpcSeedTrackMap::Reset()
{
  std::cout << "TpcSeedTrackMap: Reset() not implemented by daughter class" << std::endl;
  gSystem->Exit(1);
}
