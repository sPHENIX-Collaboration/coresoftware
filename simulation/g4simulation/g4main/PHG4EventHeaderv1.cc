#include "PHG4EventHeaderv1.h"

#include <cmath>

using namespace std;

PHG4EventHeaderv1::PHG4EventHeaderv1():
  evtseq(-9999),
  bimp(NAN),
  rplane(NAN)
{}

int
PHG4EventHeaderv1::isValid() const
{
  if (evtseq > 0)
    {
      return 1;
    }
  return 0;
}

void
PHG4EventHeaderv1::Reset()
{
  evtseq = -9999;
  bimp = NAN;
  rplane=NAN;
}

void
PHG4EventHeaderv1::identify(std::ostream& os) const
{
  os << "identify yourself: PHG4EventHeaderv1: evtseq: "
     << evtseq << ", bimp: " << bimp << ", rplane angle: " << rplane
     << endl;
}
