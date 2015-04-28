#include "PHG4CylinderGeom.h"
#include <cmath>

ClassImp(PHG4CylinderGeom)

using namespace std;

void
PHG4CylinderGeom::identify(std::ostream& os) const
{
  os << "virtual PHG4CylinderGeom"
     << endl;
  return;
}
