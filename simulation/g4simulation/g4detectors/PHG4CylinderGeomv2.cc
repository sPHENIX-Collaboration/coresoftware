#include "PHG4CylinderGeomv2.h"

#include <phparameter/PHParameters.h>

void PHG4CylinderGeomv2::identify(std::ostream& os) const
{
  os << "PHG4CylinderGeomv2: layer: " << layer
     << ", radius: " << radius
     << ", thickness: " << thickness
     << ", zmin: " << zmin
     << ", zmax: " << zmax
     << ", num scint: " << nscint
     << std::endl;
  return;
}

void PHG4CylinderGeomv2::ImportParameters(const PHParameters& param)
{
  PHG4CylinderGeomv1::ImportParameters(param);

  if (param.exist_int_param("nscint")) nscint = param.get_int_param("nscint");

  return;
}
