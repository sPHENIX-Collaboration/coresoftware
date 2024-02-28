#include "PHG4CylinderGeomv3.h"

#include <cmath>

void PHG4CylinderGeomv3::identify(std::ostream& os) const
{
  os << "PHG4CylinderGeomv3: layer: " << layer
     << ", radius: " << radius
     << ", thickness: " << thickness
     << ", zmin: " << zmin
     << ", zmax: " << zmax
     << ", num scint: " << nscint
     << ", tilt: " << tiltangle / M_PI * 180. << " deg"
     << ", phi at slat #0: " << phi_slat_zero / M_PI * 180.
     << std::endl;
  return;
}
