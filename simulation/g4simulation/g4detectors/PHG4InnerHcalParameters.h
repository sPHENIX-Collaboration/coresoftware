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

  void SetActive(const int n) {active = n;}
  int IsActive() const {return active;}

  void SetAbsorberactive(const int n) {absorberactive = n;}
  int IsAbsorberactive() const {return absorberactive;}

  void AddAbsorberHitsToTruth(const int n) {absorbertruth = n;}
  int AddAbsorberHitsToTruth() const {return absorbertruth;}

  void BlackHole(const int n) {blackhole = n;}
  int IsBlackHole() const {return blackhole;}

  void set_ncross(const int n) {ncross = n;}
  int get_ncross() const {return ncross;}

  void set_n_scinti_plates(const int n) {n_scinti_plates = n;}
  int get_n_scinti_plates() const {return n_scinti_plates;}

  void set_n_scinti_tiles(const int n) {n_scinti_tiles = n;}
  int get_n_scinti_tiles() const {return n_scinti_tiles;}

  void set_light_scint_model(const int n) {light_scint_model = n;}
  int get_light_scint_model() const {return light_scint_model;}

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

  void set_steplimits(const double x) {steplimits = x;}
  double get_steplimits() const;

  void SetLightCorrection(const double inner_radius, const double inner_corr,
			const double outer_radius, const double outer_corr);

  double get_light_balance_inner_radius() const;
  double get_light_balance_outer_radius() const;
  double get_light_balance_inner_corr() const {return light_balance_inner_corr;}
  double get_light_balance_outer_corr() const {return light_balance_outer_corr;}
  void set_place_z(const double x) {place_z = x;}
  void set_place(const double x, const double y, const double z) 
  {
    place_x = x;
    place_y = y;
    place_z = z;
  }
  double get_place_x() const;
  double get_place_y() const;
  double get_place_z() const;

  void set_rot_x(const double x) {rot_x = x;}
  void set_rot_y(const double x) {rot_y = x;}
  void set_rot_z(const double x) {rot_z = x;}
  double get_rot_x() const;
  double get_rot_y() const;
  double get_rot_z() const;

 protected:
  std::string material;
  int active;
  int absorberactive;
  int absorbertruth;
  int blackhole;
  int ncross;
  int n_scinti_plates;
  int n_scinti_tiles;
  int light_scint_model;
  double inner_radius;
  double outer_radius;
  double size_z;
  double scinti_gap;
  double tilt_angle;
  double scinti_tile_thickness;
  double scinti_gap_neighbor;
  double scinti_eta_coverage;
  double steplimits;
  double light_balance_inner_radius;
  double light_balance_inner_corr;
  double light_balance_outer_radius;
  double light_balance_outer_corr;
  double place_x;
  double place_y;
  double place_z;
  double rot_x;
  double rot_y;
  double rot_z;
};

#endif
