#include "PHG4EventHeaderv1.h"

int PHG4EventHeaderv1::isValid() const
{
  if (evtseq > 0)
  {
    return 1;
  }
  return 0;
}

void PHG4EventHeaderv1::Reset()
{
  evtseq = -9999;
  bimp = std::numeric_limits<float>::quiet_NaN();
  rplane = std::numeric_limits<float>::quiet_NaN();
}

void PHG4EventHeaderv1::identify(std::ostream& os) const
{
  os << "identify yourself: PHG4EventHeaderv1: evtseq: "
     << evtseq << ", bimp: " << bimp << ", rplane angle: " << rplane
     << std::endl;
}
