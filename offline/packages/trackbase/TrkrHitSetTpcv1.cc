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
  // not resetting properties as the hitset are used in the next events
  //  m_hitSetKey = TrkrDefs::HITSETKEYMAX;
  //  m_StartingBCO = 0;
  //  m_padIndexStart = 0;
  //  m_tBinIndexStart = 0;

  // ADC 2D array is zeroed but intentionally NOT deleted,
  // which potentially save memory allocation-deallocation ops in the next event
  for (auto& pad : m_timeFrameADCData)
  {
    std::fill(pad.begin(), pad.end(), 0);
  }
}
void TrkrHitSetTpcv1::Resize()
{
  const uint16_t n_pad = getNPads();
  const uint16_t n_tbin = getNTBins();

  if (n_pad != m_timeFrameADCData.size())
    m_timeFrameADCData.resize(n_pad);

  for (auto& pad : m_timeFrameADCData)
  {
    if (n_tbin != pad.size())
    {
      pad.reserve(n_tbin);
      pad.resize(n_tbin, 0);
    }
  }
}

void TrkrHitSetTpcv1::identify(std::ostream& os) const
{
  const uint16_t layer = TrkrDefs::getLayer(m_hitSetKey);
  const uint16_t trkrid = TrkrDefs::getTrkrId(m_hitSetKey);

  os
      << "TrkrHitSetTpcv1: "
      << "       hitsetkey " << getHitSetKey()
      << " TrkrId: " << trkrid
      << " layer: " << layer
      << " m_nPads: " << m_nPads
      << " n_tBins: " << m_nTBins
      << " m_padIndexStart: " << m_padIndexStart
      << " m_tBinIndexStart: " << m_tBinIndexStart
      << " m_StartingBCO: " << m_StartingBCO
      << std::endl;

  //  for( const auto& entry : m_hits )
  //  {
  //    std::cout << " hitkey " << entry.first << std::endl;
  //    (entry.second)->identify(os);
  //  }

  for (uint16_t i = 0; i < m_nPads; ++i)
  {
    if (m_timeFrameADCData[i].size() == 0)
    {
      os << "Pad " << i << " ADC vector is zero-sized" << std::endl;
      continue;
    }

    // skip empty pads
    if (*std::max_element(m_timeFrameADCData[i].begin(), m_timeFrameADCData[i].end()) == 0) continue;

    os << "Pad " << i << ":";

    int j = 0;
    for (const auto& adc : m_timeFrameADCData[i])
    {
      if (adc)
      {
        os << "\t[" << j << "]:" << adc;
      }
      ++j;
    }
    os << std::endl;
  }
}

TpcDefs::ADCDataType& TrkrHitSetTpcv1::getTpcADC(const uint16_t pad, const uint16_t tbin)
{
  assert(pad < m_nPads);
  assert(tbin < m_nTBins);
  assert(m_timeFrameADCData.size() == m_nPads);
  assert(m_timeFrameADCData[pad].size() == m_nTBins);

  return m_timeFrameADCData[pad][tbin];
}

const TpcDefs::ADCDataType& TrkrHitSetTpcv1::getTpcADC(const uint16_t pad, const uint16_t tbin) const
{
  assert(pad < m_nPads);
  assert(tbin < m_nTBins);
  assert(m_timeFrameADCData.size() == m_nPads);
  assert(m_timeFrameADCData[pad].size() == m_nTBins);

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

  exit(1);

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
