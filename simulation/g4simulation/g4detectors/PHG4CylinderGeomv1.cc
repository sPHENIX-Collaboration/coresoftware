#include "PHG4CylinderGeomv1.h"
#include <cmath>
#include "PHG4Parameters.h"

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

  if (param.exist_int_param("layer")) layer = param.get_int_param("layer");
  if (param.exist_double_param("radius")) radius = param.get_double_param("radius");
  if (param.exist_double_param("zmin")) zmin = param.get_double_param("zmin");
  if (param.exist_double_param("zmax")) zmax = param.get_double_param("zmax");
  if (param.exist_double_param("thickness")) thickness = param.get_double_param("thickness");

  return;
}
