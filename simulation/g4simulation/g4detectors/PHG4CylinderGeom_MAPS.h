#ifndef PHG4CylinderGeomMAPS_H__
#define PHG4CylinderGeomMAPS_H__

#include "PHG4CylinderGeom.h"
#include "TVector3.h"
#include <phool/phool.h>
#include <cmath>

class PHG4CylinderGeom_MAPS: public PHG4CylinderGeom
{
 public:

  PHG4CylinderGeom_MAPS(int layer, int stave_type, int in_Nstaves, double in_layer_nominal_radius, double in_phistep, double in_phitilt, double in_pixel_x, double in_pixel_y, double in_pixel_thickness);

  virtual ~PHG4CylinderGeom_MAPS() {}

  void identify(std::ostream& os = std::cout) const;
  TVector3 get_local_from_world_coords(int stave, int half_stave, int module, int chip, TVector3 world_location);
  TVector3 get_world_from_local_coords(int stave, int half_stave, int module, int chip, TVector3 sensor_local);
  int get_pixel_from_local_coords(TVector3 sensor_local);
  TVector3 get_local_coords_from_pixel(int NXZ);
  int get_pixel_X_from_pixel_number(int NXZ) ;
  int get_pixel_Y_from_pixel_number(int NXZ) ;

  double get_pixel_x() const {return pixel_x;}
  double get_pixel_y() const {return pixel_y;}
  double get_pixel_thickness() const {return pixel_thickness;}
  
  void set_layer(const int i) {layer = i;}
  int get_layer() const {return layer;}
  double get_radius() const {return layer_radius;}

  int get_ladder_phi_index(int stave, int half_stave, int chip);
  int get_ladder_z_index(int module, int chip);

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

  int NZ;
  int NX;

  // finding the center of a stave
  double layer_radius;
  double stave_phi_step;
  double stave_phi_tilt;

  // finding the sensor location
  double sensor_z_spacing;

  // for all layers
  double loc_sensor_in_chip[3];

  // inner barrel layers stave construction 
  // (stave_type == 0)
  double inner_loc_chip_in_module[9][3];
  double inner_loc_module_in_halfstave[3];
  double inner_loc_halfstave_in_stave[3];

  // middle barrel layers stave construction 
  // (stave_type == 1)
  double middle_loc_chip_in_module[14][3];
  double middle_loc_module_in_halfstave[4][3];
  double middle_loc_halfstave_in_stave[2][3];

  // outer barrel layers stave construction 
  // (stave_type == 2)
  double outer_loc_chip_in_module[14][3];
  double outer_loc_module_in_halfstave[7][3];
  double outer_loc_halfstave_in_stave[2][3];

  double pixel_x;
  double pixel_y;
  double pixel_thickness;
  
  ClassDef(PHG4CylinderGeom_MAPS,1)
};

#endif
