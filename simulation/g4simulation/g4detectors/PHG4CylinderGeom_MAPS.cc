#include "PHG4CylinderGeom_MAPS.h"
#include <cmath>

ClassImp(PHG4CylinderGeom_MAPS)

using namespace std;

PHG4CylinderGeom_MAPS::PHG4CylinderGeom_MAPS(int in_layer, int in_Nstaves, double in_layer_nominal_radius, double in_phistep, double in_phitilt):
  layer(in_layer),
  N_staves(in_Nstaves),
  N_half_staves(0),
  N_modules(0),
  N_chips(0),
  layer_radius(in_layer_nominal_radius),
  stave_phi_step(in_phistep),
  stave_phi_tilt(in_phitilt)
{
  // construct the geometry
  construct_geometry();

  return;
}

void
PHG4CylinderGeom_MAPS::construct_geometry()
{
  if(layer < 2)
    N_half_staves = N_staves;
  else
    N_half_staves = N_staves * 2;
  
  N_modules = 1;
  N_chips = 1;

  return;
}

void
PHG4CylinderGeom_MAPS::identify(std::ostream& os) const
{
  os << "PHG4CylinderGeom_MAPS: layer: " << layer 
     << ", layer_radius: " << layer_radius 
     << ", N_staves in layer: " << N_staves
     << ", N_half_staves in layer: " << N_half_staves
     << ", N_modules in layer: " << N_modules
     << ", N_chips in layer: " << N_chips
     << endl;
  return;
}

void PHG4CylinderGeom_MAPS::find_sensor_center(int stave_number, int half_stave_number, int module_number, int chip_number, double location[])
{
  double x_location = 0.0;
  double y_location = 0.0;
  double z_location = 0.0;

  location[0] = x_location;
  location[1] = y_location;
  location[2] = z_location;
}


