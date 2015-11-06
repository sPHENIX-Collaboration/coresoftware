#ifndef PHG4INNERHCALPARAMETERS_H
#define PHG4INNERHCALPARAMETERS_H

//#include <Geant4/globals.hh>
#include <string>
// contains parameters in G4 internal units
class PHG4InnerHcalParameters
{
 public:
  PHG4InnerHcalParameters();
  virtual ~PHG4InnerHcalParameters() {}
  void print() const;

  void set_inner_radius(const double x) {inner_radius = x;}
  double get_inner_radius() const;
  
  void set_outer_radius(const double x) {outer_radius = x;}
  double get_outer_radius() const;

 protected:
  double inner_radius;
  double outer_radius;
 public:
  double size_z;
  double scinti_gap;
  double tilt_angle;
  int n_scinti_plates;
  int n_scinti_tiles;
  double scinti_tile_thickness;
  double scinti_gap_neighbor;
  double scinti_eta_coverage;
  double place_in_x;
  double place_in_y;
  double place_in_z;
  double x_rot;
  double y_rot;
  double z_rot;
  int active;
  int absorberactive;
  int ncross;
  int blackhole;
  std::string material;
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
