#ifndef PHG4CylinderGeomMAPS_H__
#define PHG4CylinderGeomMAPS_H__

#include "PHG4CylinderGeomv4.h"
#include "TVector3.h"
#include <phool/phool.h>
#include <cmath>

class PHG4CylinderGeom_MAPS: public PHG4CylinderGeomv4
{
 public:

  PHG4CylinderGeom_MAPS(int layer, int stave_type, int in_Nstaves, double in_layer_nominal_radius, double in_phistep, double in_phitilt);

  virtual ~PHG4CylinderGeom_MAPS() {}

  void identify(std::ostream& os = std::cout) const;
  TVector3 get_world_from_local_coords(int stave, int half_stave, int module, int chip, TVector3 sensor_local);
  int get_pixel_from_local_coords(TVector3 sensor_local);
  TVector3 get_local_coords_from_pixel(int NXZ);

  void set_layer(const int i) {layer = i;}
  int get_layer() const {return layer;}
  double get_radius() const {return layer_radius;}

  void find_sensor_center(int stave_number, int half_stave_number, int module_number, int chip_number, double location[]);
  int get_N_staves() const {return N_staves;}
  int get_N_half_staves() const {return N_half_staves;}
  //int get_N_modules() const {return N_modules;}
  //int get_N_chips() const {return N_chips;}

protected:

  int layer;
  int stave_type;
  int N_staves;
  int N_half_staves;
  //int N_modules;
  //int N_chips;

  double Xsensor;
  double Zsensor;

  double pixel_x;
  double pixel_z;

  int NZ;
  int NX;

  // finding the center of a stave
  double layer_radius;
  double stave_phi_step;
  double stave_phi_tilt;

  // finding the sensor location
  double sensor_z_spacing;

  
  ClassDef(PHG4CylinderGeom_MAPS,1)
};

#endif
