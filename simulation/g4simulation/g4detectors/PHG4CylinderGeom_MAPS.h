#ifndef PHG4CylinderGeomMAPS_H__
#define PHG4CylinderGeomMAPS_H__

#include "PHG4CylinderGeomv4.h"

#include <phool/phool.h>
#include <cmath>

class PHG4CylinderGeom_MAPS: public PHG4CylinderGeomv4
{
 public:

  PHG4CylinderGeom_MAPS(int layer, int in_Nstaves, double in_layer_nominal_radius, double in_phistep, double in_phitilt);

  virtual ~PHG4CylinderGeom_MAPS() {}

  void identify(std::ostream& os = std::cout) const;
  void construct_geometry();
  void set_layer(const int i) {layer = i;}
  int get_layer() const {return layer;}
  double get_radius() const {return layer_radius;}

  void find_sensor_center(int stave_number, int half_stave_number, int module_number, int chip_number, double location[]);
  int get_N_staves() const {return N_staves;}
  int get_N_half_staves() const {return N_half_staves;}
  int get_N_modules() const {return N_modules;}
  int get_N_chips() const {return N_chips;}

protected:

  int layer;
  int N_staves;
  int N_half_staves;
  int N_modules;
  int N_chips;


  // finding the center of a stave
  double layer_radius;
  double stave_phi_step;
  double stave_phi_tilt;

  // finding the sensor location
  double sensor_z_spacing;

  
  ClassDef(PHG4CylinderGeom_MAPS,1)
};

#endif
