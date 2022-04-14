#include "TrkrClusterCrossingAssoc.h"

#include <TSystem.h>
#include <iostream>

void TrkrClusterCrossingAssoc::Reset()
{
  std::cout << "TrkrClusterCrossingAssoc: Reset() not implemented by daughter class" << std::endl;
  gSystem->Exit(1);
}
