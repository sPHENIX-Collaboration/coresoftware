#include "PHG4BlockGeom.h"

#include <ostream>  // for operator<<, endl, ostream, basic_ostream

using namespace std;

void
PHG4BlockGeom::identify(std::ostream& os) const
{
  os << "virtual PHG4BlockGeom"
     << endl;
  return;
}
