#ifndef MVTX_CYLINDERGEOMMVTX_H
#define MVTX_CYLINDERGEOMMVTX_H

#include <g4detectors/PHG4CylinderGeom.h>

#include <TVector3.h>

#include <iostream>

class CylinderGeom_Mvtx : public PHG4CylinderGeom
{
 public:
  CylinderGeom_Mvtx(
    int layer,
    int in_Nstaves,
    double in_layer_nominal_radius,
    double in_phistep,
    double in_phitilt,
    double in_phi0
  );

  //! default ctor to allow ROOT stream of this class. Implemented using c++11 feature of delegating constructors
  CylinderGeom_Mvtx()
    : CylinderGeom_Mvtx(
          /*int layer*/ 0,
          /*int in_Nstaves*/ 0,
          /*double in_layer_nominal_radius*/ 3,
          /*double in_phistep*/ 0,
          /*double in_phitilt*/ 0,
          /*double in_phi0*/ 0)
  {
  }

  virtual ~CylinderGeom_Mvtx() {}

  void identify(std::ostream& os = std::cout) const;
  TVector3 get_local_from_world_coords(int stave, int half_stave, int module, int chip, TVector3 world_location);
  TVector3 get_world_from_local_coords(int stave, int half_stave, int module, int chip, TVector3 sensor_local);

  TVector3 get_local_from_world_coords(int stave, int chip, TVector3 world_location)
  {
    return get_local_from_world_coords(stave, 0, 0, chip, world_location);
  }


  TVector3 get_world_from_local_coords(int stave, int chip, TVector3 sensor_local)
  {
    return get_world_from_local_coords(stave, 0, 0, chip, sensor_local);
  }

  bool get_pixel_from_local_coords(TVector3 sensor_local, int& iRow, int& iCol);
  int get_pixel_from_local_coords(TVector3 sensor_local);

  TVector3 get_local_coords_from_pixel(int NXZ);
  TVector3 get_local_coords_from_pixel(int iRow, int iCol);

  int get_pixel_X_from_pixel_number(int NXZ)
  {
    return NXZ % get_NX();
  }

  int get_pixel_Z_from_pixel_number(int NXZ)
  {
    return NXZ / get_NX();
  }

  int get_pixel_number_from_xbin_zbin(int xbin, int zbin) // obsolete
  {
      return xbin + zbin * get_NX();
  }

  double get_pixel_x() const { return pixel_x; }  // pitch
  double get_pixel_z() const { return pixel_z; }  // length
  double get_pixel_thickness() const { return pixel_thickness; }

  void set_layer(const int i) { layer = i; }
  int get_layer() const { return layer; }
  double get_radius() const { return layer_radius; }
  double get_stave_phi_tilt() const { return stave_phi_tilt; }
  double get_stave_phi_0() const { return stave_phi_0; }

  int get_ladder_phi_index(int stave, int half_stave, int chip) {return stave; }
  int get_ladder_z_index(int module, int chip) { return chip; }

  void find_sensor_center(int stave_number, int half_stave_number, int module_number, int chip_number, double location[]);

  int get_N_staves() const { return N_staves; }
  int get_N_half_staves() const { return N_half_staves; }

  int get_NZ() const;
  int get_NX() const;

 protected:
  int layer;
  int N_staves;
  int N_half_staves;

  // finding the center of a stave
  double layer_radius;
  double stave_phi_step;
  double stave_phi_tilt;
  double stave_phi_0;

  // for all layers
  double loc_sensor_in_chip[3];

  // inner barrel layers stave construction
  double inner_loc_chip_in_module[9][3];
  double inner_loc_module_in_halfstave[3];
  double inner_loc_halfstave_in_stave[3];

  double pixel_x;
  double pixel_z;
  double pixel_thickness;

  ClassDef(CylinderGeom_Mvtx, 2)
};

#endif
