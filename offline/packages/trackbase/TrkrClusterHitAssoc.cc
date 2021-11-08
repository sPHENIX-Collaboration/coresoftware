#include "TrkrClusterHitAssoc.h"

#include <TSystem.h>
#include <iostream>

void TrkrClusterHitAssoc::Reset()
{
  std::cout << "TrkrClusterHitAssoc: Reset() not implemented by daughter class" << std::endl;
  gSystem->Exit(1);
}
