#include "TrkrHitv1.h"

void TrkrHitv1::CopyFrom(const TrkrHit& source)
{
  // do nothing if copying onto oneself
  if (this == &source)
  {
    return;
  }

  // parent class method
  TrkrHit::CopyFrom(source);

  // copy adc
  setAdc(source.getAdc());
}

unsigned int TrkrHitv1::getAdc() const
{
  return m_adc;
}
