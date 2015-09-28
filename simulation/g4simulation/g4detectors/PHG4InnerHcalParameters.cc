#include "PHG4InnerHcalParameters.h"

#include <Geant4/G4SystemOfUnits.hh>

#include <iostream>

using namespace std;

PHG4InnerHcalParameters::PHG4InnerHcalParameters():
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
  steplimits(NAN)
{}

void
PHG4InnerHcalParameters::print() const
{
  cout << "Inner Radius: " << inner_radius/cm << endl;
  cout << "Outer Radius: " << outer_radius/cm << endl;
  cout << "Size Z: " << size_z/cm << endl;
  cout << "Scintillator Gap: " << scinti_gap/cm << endl;
  cout << "Tilt Angle: " << tilt_angle/deg << endl;
  cout << "Crossings: " << ncross << endl;
  return;
}
