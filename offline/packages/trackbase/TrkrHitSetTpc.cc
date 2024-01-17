/**
 * @file trackbase/TrkrHitSetTpc.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of TrkrHitSetTpc
 */
#include "TrkrHitSetTpc.h"
#include "TrkrHit.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>  // for exit
#include <iostream>
#include <type_traits>  // for __decay_and_strip<>::__type

void TrkrHitSetTpc::identify(std::ostream& os) const
{
  os
      << "TrkrHitSetTpc: "
      << "       hitsetkey " << getHitSetKey()
      << std::endl;
}

TpcDefs::ADCDataType& TrkrHitSetTpc::getTpcADC(TrkrDefs::hitkey key)
{
  const uint16_t pad = TpcDefs ::getPad(key);
  const uint16_t tbin = TpcDefs ::getTBin(key);

  return getTpcADC(pad, tbin);
}

const TpcDefs::ADCDataType& TrkrHitSetTpc::getTpcADC(TrkrDefs::hitkey key) const
{
  const uint16_t pad = TpcDefs ::getPad(key);
  const uint16_t tbin = TpcDefs ::getTBin(key);

  return getTpcADC(pad, tbin);
}
