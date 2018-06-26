/**
 * @file mvtx/MvtxHit.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of Mvtx hit object
 */
#include "MvtxHit.h"
#include "MvtxDefs.h"

MvtxHit::MvtxHit()
  : TrkrHit()
{
}

void 
MvtxHit::identify(std::ostream& os) const
{
  os << "MvtxHit with key:" << getKey() 
     << std::endl;
}

void 
MvtxHit::Reset()
{
  TrkrHit::Reset();
}

int 
MvtxHit::isValid() const
{
  // valid if the key is not equal to the default value
  return getKey() != TrkrDefs::HITKEYMAX;
}
