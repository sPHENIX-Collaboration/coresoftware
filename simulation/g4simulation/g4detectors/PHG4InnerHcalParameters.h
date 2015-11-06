#ifndef PHG4INNERHCALPARAMETERS_H
#define PHG4INNERHCALPARAMETERS_H

#include <string>

// contains parameters in our units,
// convert to G4 units inside get access methods

class PHG4InnerHcalParameters
{
 public:
  PHG4InnerHcalParameters();
  virtual ~PHG4InnerHcalParameters() {}
  void print() const;

  void set_material(const std::string mat) {material = mat;}
  std::string get_material() const {return material;}

  void set_ncross(const int n) {ncross = n;}
  int get_ncross() const {return ncross;}

  void set_n_scinti_plates(const int n) {n_scinti_plates = n;}
  int get_n_scinti_plates() const {return n_scinti_plates;}

  void set_n_scinti_tiles(const int n) {n_scinti_tiles = n;}
  int get_n_scinti_tiles() const {return n_scinti_tiles;}

  void set_inner_radius(const double x) {inner_radius = x;}
  double get_inner_radius() const;
  
  void set_outer_radius(const double x) {outer_radius = x;}
  double get_outer_radius() const;

  void set_size_z(const double x) {size_z = x;}
  double get_size_z() const;

  void set_scinti_gap(const double x) {scinti_gap = x;}
  double get_scinti_gap() const;

  void set_tilt_angle(const double x, const int resetncross = 1);
  double get_tilt_angle() const;

  void set_scinti_tile_thickness(const double x) {scinti_tile_thickness = x;}
  double get_scinti_tile_thickness() const;

  void set_scinti_gap_neighbor(const double x) {scinti_gap_neighbor = x;}
  double get_scinti_gap_neighbor() const;

  void set_scinti_eta_coverage(const double x) {scinti_eta_coverage = x;}
  double get_scinti_eta_coverage() const {return scinti_eta_coverage;}

 protected:
  std::string material;
  int ncross;
  int n_scinti_plates;
  int n_scinti_tiles;
  double inner_radius;
  double outer_radius;
  double size_z;
  double scinti_gap;
  double tilt_angle;
  double scinti_tile_thickness;
  double scinti_gap_neighbor;
  double scinti_eta_coverage;
 public:
  double place_in_x;
  double place_in_y;
  double place_in_z;
  double x_rot;
  double y_rot;
  double z_rot;
  int active;
  int absorberactive;
  int blackhole;
  double steplimits;
  bool  light_scint_model;
  bool light_balance;
  double light_balance_inner_radius;
  double light_balance_inner_corr;
  double light_balance_outer_radius;
  double light_balance_outer_corr;
  int absorbertruth;
};

#endif
