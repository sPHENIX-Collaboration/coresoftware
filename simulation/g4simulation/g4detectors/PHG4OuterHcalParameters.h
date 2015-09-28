#ifndef PHG4OUTERHCALPARAMETERS_H
#define PHG4OUTERHCALPARAMETERS_H

#include <Geant4/globals.hh>

// contains parameters in G4 internal units
class PHG4OuterHcalParameters
{
 public:
  PHG4OuterHcalParameters();
  virtual ~PHG4OuterHcalParameters() {}
  void print() const;
  G4double inner_radius;
  G4double outer_radius;
  G4double size_z;
  G4double scinti_gap;
  G4double tilt_angle;
  G4int n_scinti_plates;
  G4int n_scinti_tiles;
  G4double scinti_tile_thickness;
  G4double scinti_gap_neighbor;
  G4double scinti_eta_coverage;
  G4double place_in_x;
  G4double place_in_y;
  G4double place_in_z;
  G4double x_rot;
  G4double y_rot;
  G4double z_rot;
  G4int active;
  G4int absorberactive;
  G4int ncross;
  G4int blackhole;
  G4String material;
  G4double steplimits;
  bool enable_field_checker;
  bool  light_scint_model;
  bool light_balance;
  G4double light_balance_inner_radius;
  G4double light_balance_inner_corr;
  G4double light_balance_outer_radius;
  G4double light_balance_outer_corr;
  G4double magnet_cutout;
  G4int magnet_cutout_first_scinti;
};

#endif
