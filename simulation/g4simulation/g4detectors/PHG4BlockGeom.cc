#include "PHG4BlockGeom.h"
#include <cmath>

ClassImp(PHG4BlockGeom)

using namespace std;

void
PHG4BlockGeom::identify(std::ostream& os) const
{
  os << "virtual PHG4BlockGeom"
     << endl;
  return;
}
