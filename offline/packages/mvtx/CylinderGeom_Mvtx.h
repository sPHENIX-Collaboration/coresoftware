#ifndef MVTX_CYLINDERGEOMMVTX_H
#define MVTX_CYLINDERGEOMMVTX_H

#include <g4detectors/PHG4CylinderGeom.h>

#include <TVector3.h>

#include <iostream>

class CylinderGeom_Mvtx : public PHG4CylinderGeom
{
 public:
  CylinderGeom_Mvtx(int layer, int stave_type, int in_Nstaves, double in_layer_nominal_radius, double in_phistep, double in_phitilt, double in_phi0, double in_pixel_x, double in_pixel_z, double in_pixel_thickness);

  //! default ctor to allow ROOT stream of this class. Implemented using c++11 feature of delegating constructors
  CylinderGeom_Mvtx()
    : CylinderGeom_Mvtx(
          /*int layer*/ 0,
          /*int stave_type*/ 0,
          /*int in_Nstaves*/ 0,
          /*double in_layer_nominal_radius*/ 3,
          /*double in_phistep*/ 0,
          /*double in_phitilt*/ 0,
          /*double in_phi0*/ 0,
          /*double in_pixel_x*/ 20e-4,
          /*double in_pixel_z*/ 20e-4,
          /*double in_pixel_thickness*/ 18e-4)
  {
  }

  virtual ~CylinderGeom_Mvtx() {}

  void identify(std::ostream& os = std::cout) const;
  TVector3 get_local_from_world_coords(int stave, int half_stave, int module, int chip, TVector3 world_location);
  TVector3 get_world_from_local_coords(int stave, int half_stave, int module, int chip, TVector3 sensor_local);
  int get_pixel_from_local_coords(TVector3 sensor_local);
  TVector3 get_local_coords_from_pixel(int NXZ);
  int get_pixel_X_from_pixel_number(int NXZ);
  int get_pixel_Z_from_pixel_number(int NXZ);
  int get_pixel_number_from_xbin_zbin(int xbin, int zbin);

  double get_pixel_x() const { return pixel_x; }  // pitch
  double get_pixel_z() const { return pixel_z; }  // length
  double get_pixel_thickness() const { return pixel_thickness; }

  void set_layer(const int i) { layer = i; }
  int get_layer() const { return layer; }
  double get_radius() const { return layer_radius; }
  double get_stave_phi_tilt() const { return stave_phi_tilt; }

  int get_ladder_phi_index(int stave, int half_stave, int chip);
  int get_ladder_z_index(int module, int chip);

  void find_sensor_center(int stave_number, int half_stave_number, int module_number, int chip_number, double location[]);
  int get_N_staves() const { return N_staves; }
  int get_N_half_staves() const { return N_half_staves; }
  //int get_N_modules() const {return N_modules;}
  //int get_N_chips() const {return N_chips;}

  int get_NZ() const { return (int) (Zsensor / (pixel_z)); }
  int get_NX() const { return (int) (Xsensor / (pixel_x)); }

 protected:
  int layer;
  int stave_type;
  int N_staves;
  int N_half_staves;

  double Xsensor;
  double Zsensor;

  // finding the center of a stave
  double layer_radius;
  double stave_phi_step;
  double stave_phi_tilt;
  double stave_phi_0;

  // finding the sensor location
  //  double sensor_z_spacing;

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
  double pixel_z;
  double pixel_thickness;

  ClassDef(CylinderGeom_Mvtx, 1)
};

#endif
