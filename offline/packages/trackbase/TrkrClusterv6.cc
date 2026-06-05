/**
 * @file trackbase/TrkrClusterv6.cc
 * @author Ishan Goel
 * @date May 2026
 * @brief Implementation of TrkrClusterv6
 */
#include "TrkrClusterv6.h"

#include <cmath>
#include <utility>  // for swap

namespace
{
  // square convenience function
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }
}  // namespace

TrkrClusterv6::TrkrClusterv6()
  : m_subsurfkey(TrkrDefs::SUBSURFKEYMAX)
  , m_phierr(0)
  , m_zerr(0)
  , m_adc(0)
  , m_maxadc(0)
  , m_cenadc(0)
  , m_padcen(0)
  , m_tbincen(0)
  , m_padmax(0)
  , m_tbinmax(0)
  , m_rsize(0)
  , m_phisize(0)
  , m_zsize(0)
  , m_overlap(0)
  , m_edge(0)
  , m_sledge(0)
  , m_sredge(0)
  , m_tledge(0)
  , m_tredge(0)
  , m_dledge(0)
  , m_dredge(0)
  , m_hledge(0)
  , m_hredge(0)
  , m_slmix(0)
  , m_srmix(0)
  , m_tlmix(0)
  , m_trmix(0)
  , m_phibinlo(0)
  , m_phibinhi(0)
  , m_tbinlo(0)
  , m_tbinhi(0)
  , m_padphase(0)
  , m_tbinphase(0)
{
  for (float& i : m_local)
  {
    i = NAN;
  }
}

void TrkrClusterv6::identify(std::ostream& os) const
{
  os << "---TrkrClusterv6--------------------" << std::endl;

  os << " (rphi,z) =  (" << getLocalX();
  os << ", " << getLocalY() << ") cm ";

  os << " valid = " << isValid() << std::endl;

  os << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int TrkrClusterv6::isValid() const
{
  for (int i = 0; i < 2; ++i)
  {
    if (std::isnan(getPosition(i)))
    {
      return 0;
    }
  }
  if (m_adc == 0xFFFF)
  {
    return 0;
  }

  return 1;
}

void TrkrClusterv6::CopyFrom(const TrkrCluster& source)
{
  // do nothing if copying onto oneself
  if (this == &source)
  {
    return;
  }

  // parent class method
  TrkrCluster::CopyFrom(source);

  setLocalX(source.getLocalX());
  setLocalY(source.getLocalY());
  setSubSurfKey(source.getSubSurfKey());
  setAdc(source.getAdc());
  setMaxAdc(source.getMaxAdc());
  setCenAdc(source.getCenAdc());
  setPadCen(source.getPadCen());
  setTBinCen(source.getTBinCen());
  setPadMax(source.getPadMax());
  setTBinMax(source.getTBinMax());
  setPhiError(source.getRPhiError());
  setZError(source.getZError());
  setRSize(source.getRSize());
  setPhiSize(source.getPhiSize());
  setZSize(source.getZSize());
  setOverlap(source.getOverlap());
  setEdge(source.getEdge());
  setSLEdge(source.getSLEdge());
  setSREdge(source.getSREdge());
  setTLEdge(source.getTLEdge());
  setTREdge(source.getTREdge());
  setDLEdge(source.getDLEdge());
  setDREdge(source.getDREdge());
  setHLEdge(source.getHLEdge());
  setHREdge(source.getHREdge());
  setSLMix(source.getSLMix());
  setSRMix(source.getSRMix());
  setTLMix(source.getTLMix());
  setTRMix(source.getTRMix());
  setPhiBinLo(source.getPhiBinLo());
  setPhiBinHi(source.getPhiBinHi());
  setTBinLo(source.getTBinLo());
  setTBinHi(source.getTBinHi());
  setPadPhase(source.getPadPhase());
  setTBinPhase(source.getTBinPhase());
}
