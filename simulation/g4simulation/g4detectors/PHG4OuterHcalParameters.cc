#include "PHG4OuterHcalParameters.h"

#include <Geant4/G4SystemOfUnits.hh>

PHG4OuterHcalParameters::PHG4OuterHcalParameters():
  inner_radius(116 * cm),
  outer_radius(136 * cm),
  size_z(175.94 * 2 * cm),
  scinti_gap(0.85 * cm),
  tilt_angle(32.5 * deg),
  n_scinti_plates(5 * 64),
  n_scinti_tiles(12),
  scinti_tile_thickness(0.7*cm),
  scinti_gap_neighbor(0.1*cm),
  scinti_eta_coverage(1.1),
  place_in_x(0 * cm),
  place_in_y(0 * cm),
  place_in_z(0 * cm),
  x_rot(0*deg),
  y_rot(0*deg),
  z_rot(0*deg),
  active(0),
  absorberactive(0),
  ncross(0),
  blackhole(0),
  material("SS310"),
  steplimits(NAN),
  enable_field_checker(false),
  light_scint_model(true),
  light_balance(false),
  light_balance_inner_radius(0.0),
  light_balance_inner_corr(1.0),
  light_balance_outer_radius(10.0),
  light_balance_outer_corr(1.0)
{}
