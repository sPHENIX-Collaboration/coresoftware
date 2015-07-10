#include "PHG4InnerHcalParameters.h"

#include <Geant4/G4SystemOfUnits.hh>

PHG4InnerHcalParameters::PHG4InnerHcalParameters():
  inner_radius(116 * cm),
  outer_radius(136 * cm),
  size_z(175.94 * 2 * cm),
  scinti_gap(0.85 * cm),
  tilt_angle(32.5 * deg),
  n_scinti_plates(5 * 64),
  n_scinti_tiles(11),
  scinti_tile_thickness(0.7*cm),
  scinti_gap_neighbor(0.2*cm),
  scinti_eta_coverage(1.1),
  place_in_x(0 * cm),
  place_in_y(0 * cm),
  place_in_z(0 * cm),
  x_rot(0*deg),
  y_rot(0*deg),
  z_rot(0*deg),
  active(0),
  absorberactive(0),
  blackhole(0),
  material("SS310")
{}
