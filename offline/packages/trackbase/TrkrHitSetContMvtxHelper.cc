/**
 * @file trackbase/TrkrHitSetContMvtxHelper.cc
 * @author Yasser Corrales Morales <ycmorales@bnl.gov>
 * @date Febraury 2025
 * base class for Mvtx hitsetkey container per strobe
 */

#include "TrkrHitSetContMvtxHelper.h"

#include <Rtypes.h>
#include <TSystem.h>

#include <cstdlib>

void TrkrHitSetContMvtxHelper::Reset()
{
  std::cout << "TrkrHitSetContMvtxHelper: Reset() not implemented by daughter class" << std::endl;
  gSystem->Exit(1);
}
