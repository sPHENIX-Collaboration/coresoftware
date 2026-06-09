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
  constexpr T square(const T& x)
  {
    return x * x;
  }
}  // namespace

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
