#include "PHG4CylinderGeomv2.h"

ClassImp(PHG4CylinderGeomv2)

using namespace std;

PHG4CylinderGeomv2::PHG4CylinderGeomv2():
  nscint(-9999)
{
  return;
}

void
PHG4CylinderGeomv2::identify(std::ostream& os) const
{
  os << "PHG4CylinderGeomv2: layer: " << layer 
     << ", radius: " << radius 
     << ", thickness: " << thickness
     << ", zmin: " << zmin 
     << ", zmax: " << zmax 
     << ", num scint: " << nscint
     << endl;
  return;
}
