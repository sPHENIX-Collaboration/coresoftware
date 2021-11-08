#include "TrkrClusterIterationMap.h"

#include <TSystem.h>
#include <iostream>

void TrkrClusterIterationMap::Reset()
{
  std::cout << "TrkrClusterIterationMap: Reset() not implemented by daughter class" << std::endl;
  gSystem->Exit(1);
}
