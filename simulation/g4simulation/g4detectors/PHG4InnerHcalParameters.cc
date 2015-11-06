#include "PHG4InnerHcalParameters.h"

#include <Geant4/G4SystemOfUnits.hh>

#include <cmath>
#include <iostream>

using namespace std;

PHG4InnerHcalParameters::PHG4InnerHcalParameters():
  material("SS310"),
  ncross(4),
  n_scinti_plates(5 * 64),
  n_scinti_tiles(12),
  inner_radius(116),
  outer_radius(136),
  size_z(175.94 * 2),
  scinti_gap(0.85),
  tilt_angle(32.5),
  scinti_tile_thickness(0.7),
  scinti_gap_neighbor(0.1),
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
  steplimits(NAN),
  light_scint_model(true),
  light_balance(false),
  light_balance_inner_radius(0.0),
  light_balance_inner_corr(1.0),
  light_balance_outer_radius(10.0),
  light_balance_outer_corr(1.0),
  absorbertruth(0)
{}

void
PHG4InnerHcalParameters::print() const
{
  cout << "Material: " << material << endl;
  cout << "Crossings: " << ncross << endl;
  cout << "Plates: " << n_scinti_plates << endl;
  cout << "Tiles: " << n_scinti_tiles << endl;
  cout << "Inner Radius: " << inner_radius << endl;
  cout << "Outer Radius: " << outer_radius << endl;
  cout << "Size Z: " << size_z << endl;
  cout << "Tilt Angle: " << tilt_angle << endl;
  cout << "Scintillator Thickness: " << scinti_tile_thickness << endl;
  cout << "Scintillator Gap: " << scinti_gap << endl;
  cout << "Scintillator Gap to neighbor: " << scinti_gap_neighbor << endl;
  cout << "Eta coverage: " << scinti_eta_coverage << endl;
  cout << "Light Balance Inner Radius: " << light_balance_inner_radius << endl;
  cout << "Light Balance Inner Corr: " << light_balance_inner_corr << endl;
  cout << "Light Balance Outer Radius: " << light_balance_outer_radius << endl;
  cout << "Light Balance Outer Corr: " << light_balance_outer_corr << endl;
  return;
}

double
PHG4InnerHcalParameters::get_inner_radius() const 
{
  return inner_radius*cm;
}
double
PHG4InnerHcalParameters::get_outer_radius() const 
{
  return outer_radius*cm;
}

double
PHG4InnerHcalParameters::get_size_z() const
{
  return size_z*cm;
}

double
PHG4InnerHcalParameters::get_scinti_gap() const
{
return scinti_gap*cm;
}

void
PHG4InnerHcalParameters::set_tilt_angle(const double x, const int resetncross) 
{
  tilt_angle = x;
  if (resetncross)
    {
      ncross = 0;
    }
}

double
PHG4InnerHcalParameters::get_tilt_angle() const
{
  return tilt_angle*deg;
}

double
PHG4InnerHcalParameters::get_scinti_tile_thickness() const
{
  return scinti_tile_thickness*cm;
}

double
PHG4InnerHcalParameters::get_scinti_gap_neighbor() const
{
  return scinti_gap_neighbor*cm;
}

