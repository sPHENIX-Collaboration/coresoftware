/**
 * @file trackbase/TrkrHitSetTpcv1.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of TrkrHitSetTpcv1
 */
#include "TrkrHitSetTpcv1.h"
#include "TrkrHit.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>  // for exit
#include <iostream>
#include <type_traits>  // for __decay_and_strip<>::__type

void TrkrHitSetTpcv1::Reset()
{
  m_hitSetKey = TrkrDefs::HITSETKEYMAX;

  for (auto& pad : m_timeFrameADCData)
  {
    std::fill(pad.begin(), pad.end(), 0);
  }
}
void TrkrHitSetTpcv1::Resize(const unsigned int n_pad, const unsigned int n_tbin)
{
  m_timeFrameADCData.resize(n_pad);

  for (auto& pad : m_timeFrameADCData)
  {
    pad.resize(n_tbin);
  }

  Reset();
}

void TrkrHitSetTpcv1::identify(std::ostream& os) const
{
  const unsigned int layer = TrkrDefs::getLayer(m_hitSetKey);
  const unsigned int trkrid = TrkrDefs::getTrkrId(m_hitSetKey);

  os
      << "TrkrHitSetTpcv1: "
      << "       hitsetkey " << getHitSetKey()
      << " TrkrId " << trkrid
      << " layer " << layer
      << " m_nPads: " << m_nPads
      << " n_tBins: " << n_tBins
      << std::endl;

  //  for( const auto& entry : m_hits )
  //  {
  //    std::cout << " hitkey " << entry.first << std::endl;
  //    (entry.second)->identify(os);
  //  }

  for (unsigned int i = 0; i < m_nPads; ++i)
  {
    os << "Pad " << i << ":";

    for (const auto& adc : m_timeFrameADCData[i])
    {
      os << "\t" << adc;
    }
    os << std::endl;
  }
}

TpcDefs::ADCDataType& TrkrHitSetTpcv1::getTpcADC(const uint16_t pad, const uint16_t tbin)
{
  assert(pad < m_nPads);
  assert(tbin < n_tBins);

  return m_timeFrameADCData[pad][tbin];
}

const TpcDefs::ADCDataType& TrkrHitSetTpcv1::getTpcADC(const uint16_t pad, const uint16_t tbin) const
{
  assert(pad < m_nPads);
  assert(tbin < n_tBins);

  return m_timeFrameADCData[pad][tbin];
}

void TrkrHitSetTpcv1::removeHit(TrkrDefs::hitkey key)
{
  getTpcADC(key) = 0;
}

TrkrHitSetTpcv1::ConstIterator
TrkrHitSetTpcv1::addHitSpecificKey(const TrkrDefs::hitkey key, TrkrHit* hit)
{
  std::cout << __PRETTY_FUNCTION__
            << " : This function is deprecated! Please use getTpcADC(TrkrDefs::hitkey key)" << std::endl;

  if (hit)
  {
    getTpcADC(key) = hit->getAdc();
    delete hit;
  }

  return TrkrHitSetTpc::addHitSpecificKey(key, hit);
}

TrkrHit*
TrkrHitSetTpcv1::getHit(const TrkrDefs::hitkey key) const
{
  std::cout << __PRETTY_FUNCTION__
            << " : This function is deprecated! Please use getTpcADC(TrkrDefs::hitkey key)" << std::endl;

  exit(1);

  return TrkrHitSetTpc::getHit(key);
}

TrkrHitSetTpcv1::ConstRange
TrkrHitSetTpcv1::getHits() const
{
  std::cout << __PRETTY_FUNCTION__
            << " : This function is deprecated! Please use getTpcADC(TrkrDefs::hitkey key)" << std::endl;

  exit(1);

  return TrkrHitSetTpc::getHits();
}
