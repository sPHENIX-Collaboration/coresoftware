#include "PHG4CylinderGeom.h"

#include <ostream>  // for operator<<, endl, ostream, basic_ostream

using namespace std;

void
PHG4CylinderGeom::identify(std::ostream& os) const
{
  os << "virtual PHG4CylinderGeom"
     << endl;
  return;
}
