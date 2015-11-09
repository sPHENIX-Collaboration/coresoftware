#include "PHG4InnerHcalParameters.h"

#include <Geant4/G4SystemOfUnits.hh>

#include <cmath>
#include <iostream>

using namespace std;

PHG4InnerHcalParameters::PHG4InnerHcalParameters():
  material("SS310"),
  active(0),
  absorberactive(0),
  absorbertruth(0),
  blackhole(0),
  ncross(4),
  n_scinti_plates(5 * 64),
  n_scinti_tiles(12),
  light_scint_model(1),
  inner_radius(116),
  outer_radius(136),
  size_z(175.94 * 2),
  scinti_gap(0.85),
  tilt_angle(32.5),
  scinti_tile_thickness(0.7),
  scinti_gap_neighbor(0.1),
  scinti_eta_coverage(1.1),
  steplimits(NAN),
  light_balance_inner_radius(NAN),
  light_balance_inner_corr(NAN),
  light_balance_outer_radius(NAN),
  light_balance_outer_corr(NAN),
  place_x(0.),
  place_y(0.),
  place_z(0.),
  rot_x(0.),
  rot_y(0.),
  rot_z(0.)
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
  cout << "Max Step limit: " << steplimits << endl;
  cout << "Eta coverage: " << scinti_eta_coverage << endl;
  cout << "Light Balance Inner Radius: " << light_balance_inner_radius << endl;
  cout << "Light Balance Inner Corr: " << light_balance_inner_corr << endl;
  cout << "Light Balance Outer Radius: " << light_balance_outer_radius << endl;
  cout << "Light Balance Outer Corr: " << light_balance_outer_corr << endl;
  cout << "Placement in X: " << place_x << endl;
  cout << "Placement in Y: " << place_y << endl;
  cout << "Placement in Z: " << place_z << endl;
  cout << "Rotation in X: " << rot_x << endl;
  cout << "Rotation in Y: " << rot_y << endl;
  cout << "Rotation in Z: " << rot_z << endl;
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

double
PHG4InnerHcalParameters::get_steplimits() const
{
  return steplimits*cm;
}

void
PHG4InnerHcalParameters::SetLightCorrection(const double inner_radius, const double inner_corr,const double outer_radius, const double outer_corr)
{
  light_balance_inner_radius = inner_radius;
  light_balance_inner_corr = inner_corr;
  light_balance_outer_radius = outer_radius;
  light_balance_outer_corr = outer_corr;
}

double 
PHG4InnerHcalParameters::get_light_balance_inner_radius() const
{
  return light_balance_inner_radius*cm;
}

double 
PHG4InnerHcalParameters::get_light_balance_outer_radius() const
{
  return light_balance_outer_radius*cm;
}


double
PHG4InnerHcalParameters::get_place_x() const
{
  return place_x*cm;
}

double
PHG4InnerHcalParameters::get_place_y() const
{
  return place_y*cm;
}

double
PHG4InnerHcalParameters::get_place_z() const
{
  return place_z*cm;
}

double
PHG4InnerHcalParameters::get_rot_x() const
{
  return rot_x*deg;
}

double
PHG4InnerHcalParameters::get_rot_y() const
{
  return rot_y*deg;
}

double
PHG4InnerHcalParameters::get_rot_z() const
{
  return rot_z*deg;
}
