#include "PHG4VtxPoint.h"


ClassImp(PHG4VtxPoint)

using namespace std;

void
PHG4VtxPoint::identify(ostream& os) const
{
  os << "virtual PHG4VtxPoint base class"
     << endl;
}

bool
PHG4VtxPoint::operator== (const PHG4VtxPoint& p) const
{
  if (p.get_x() == get_x() && p.get_y() == get_y() &&
      p.get_z() == get_z() && p.get_t() == get_t())
    {
      return true;
    }
  return false;
}
 
