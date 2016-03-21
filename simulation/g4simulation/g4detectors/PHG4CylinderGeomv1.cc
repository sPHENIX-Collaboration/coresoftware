#include "PHG4CylinderGeomv1.h"
#include <cmath>
#include "PHG4Parameters.h"

ClassImp(PHG4CylinderGeomv1)

using namespace std;

PHG4CylinderGeomv1::PHG4CylinderGeomv1():
  layer(-1),
  radius(NAN),
  zmin(NAN),
  zmax(NAN),
  thickness(NAN)
{
  return;
}

void
PHG4CylinderGeomv1::identify(std::ostream& os) const
{
  os << "PHG4CylinderGeomv1: layer: " << layer 
     << ", radius: " << radius 
     << ", thickness: " << thickness
     << ", zmin: " << zmin 
     << ", zmax: " << zmax 
     << endl;
  return;
}


void
PHG4CylinderGeomv1::ImportParameters(const PHG4Parameters & param)
{
  PHG4CylinderGeom::ImportParameters(param);

  return;
}
