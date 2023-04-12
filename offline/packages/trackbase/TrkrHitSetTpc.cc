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
  std::pair<uint16_t, uint16_t> local_phi_t = getLocalPhiTBin(key);

  return getTpcADC(local_phi_t.first, local_phi_t.second);
}

const TpcDefs::ADCDataType& TrkrHitSetTpc::getTpcADC(TrkrDefs::hitkey key) const
{
  std::pair<uint16_t, uint16_t> local_phi_t = getLocalPhiTBin(key);

  return getTpcADC(local_phi_t.first, local_phi_t.second);
}

std::pair<uint16_t, uint16_t> TrkrHitSetTpc::getLocalPhiTBin(TrkrDefs::hitkey key) const
{
  const uint16_t pad = TpcDefs ::getPad(key);
  const uint16_t tbin = TpcDefs ::getTBin(key);
  const uint8_t side = TpcDefs::getSide(getHitSetKey());

  if (side == 0)
  {
    const uint16_t local_pad = pad - getPadIndexStart();
    const uint16_t local_tbin = tbin - getTBinIndexStart();

    assert(local_pad < getNPads());
    assert(local_tbin < getNTBins());

    return std::make_pair(local_pad, local_tbin);
  }
  else
  {
    const uint16_t local_pad = getNPads() - 1 - pad + getPadIndexStart();
    const uint16_t local_tbin = -tbin + getTBinIndexStart();

    assert(local_pad < getNPads());
    assert(local_tbin < getNTBins());

    return std::make_pair(local_pad, local_tbin);
  }
}

TrkrDefs::hitkey TrkrHitSetTpc::getHitKeyfromLocalBin(
    const uint16_t local_pad,
    const uint16_t local_tbin)
    const
{
  const uint8_t side = TpcDefs::getSide(getHitSetKey());

  if (side == 0)
  {
    const uint16_t pad = local_pad + getPadIndexStart();
    const uint16_t tbin = local_tbin + getTBinIndexStart();

    return TpcDefs::genHitKey(pad, tbin);
  }
  else
  {
    const uint16_t pad = getNPads() - 1 - local_pad + getPadIndexStart();
    const uint16_t tbin = -local_tbin + getTBinIndexStart();

    return TpcDefs::genHitKey(pad, tbin);
  }
}
