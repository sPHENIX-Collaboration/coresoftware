#include "TrkrHitv3.h"

void TrkrHitv3::CopyFrom(const TrkrHit& source)
{
  // do nothing if copying onto oneself
  if (this == &source)
  {
    return;
  }

  // parent class method
  TrkrHitv2::CopyFrom(source);

  // copy timing information
  setFPHXBCO(source.getFPHXBCO());
  setBCO(source.getBCO());
}
