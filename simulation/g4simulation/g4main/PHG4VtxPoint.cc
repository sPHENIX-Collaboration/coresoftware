#include "PHG4VtxPoint.h"

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
  // note that ID is not compared here, per algorithm requirement in PHG4TruthInfoContainer::AddPrimaryVertex

  if (p.get_x() == get_x() && p.get_y() == get_y() &&
      p.get_z() == get_z() && p.get_t() == get_t())
    {
      return true;
    }
  return false;
}
 
