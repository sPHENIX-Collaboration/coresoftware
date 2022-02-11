// $$Id: PHG4CylinderGeom_Spacalv2.cc,v 1.3 2014/08/28 22:18:35 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.3 $$
 * \date $$Date: 2014/08/28 22:18:35 $$
 */

#include "PHG4CylinderGeom_Spacalv2.h"

#include <phparameter/PHParameters.h>

#include <Geant4/G4PhysicalConstants.hh>

#include <CLHEP/Units/SystemOfUnits.h>  // for twopi, halfpi, pi

#include <cmath>
#include <cstdlib>  // for exit
#include <iostream>

using namespace std;

PHG4CylinderGeom_Spacalv2::PHG4CylinderGeom_Spacalv2()
{
  SetDefault();
}

void PHG4CylinderGeom_Spacalv2::identify(std::ostream& os) const
{
  os << "PHG4CylinderGeom_Spacalv2: layer: " << layer  //
     << ", radius: " << radius                         //
     << ", thickness: " << thickness                   //
     << ", zmin: " << zmin                             //
     << ", zmax: " << zmax <<                          //
      ", num scint: " << nscint

     << endl;
  return;
}

void PHG4CylinderGeom_Spacalv2::Print(Option_t* opt) const
{
  PHG4CylinderGeom_Spacalv1::Print(opt);

  cout << "\t"
       << "is_azimuthal_seg_visible() = " << is_azimuthal_seg_visible()
       << endl;
  cout << "\t"
       << "azimuthal_tilt() = " << get_azimuthal_tilt() << endl;
  cout << "\t"
       << "get_polar_taper_ratio() = " << get_polar_taper_ratio()
       << endl;
  cout << "\t"
       << "get_sec_azimuthal_width() = " << get_sec_azimuthal_width()
       << endl;
  cout << "\t"
       << "get_sec_depth() = " << get_sec_depth() << endl;
  cout << "\t"
       << "get_block_width() = " << get_block_width() << endl;
  cout << "\t"
       << "get_block_depth() = " << get_block_depth() << endl;
  cout << "\t"
       << "get_assembly_spacing() = " << get_assembly_spacing() << endl;
  cout << "\t"
       << "get_reg_fiber_grid_distance_taper() = "
       << get_reg_fiber_grid_distance_taper() << " = sqrt(3)*"
       << get_reg_fiber_grid_distance_taper() / sqrt(3.) << endl;
  cout << "\t"
       << "get_reg_fiber_grid_distance_nontaper() = "
       << get_reg_fiber_grid_distance_nontaper() << endl;
}

void PHG4CylinderGeom_Spacalv2::SetDefault()
{
  PHG4CylinderGeom_Spacalv1::SetDefault();

  azimuthal_n_sec = 256;
  azimuthal_tilt = 0;
  azimuthal_seg_visible = false;
  polar_taper_ratio = 1 + 1.1 / 42.;
  assembly_spacing = 0.0001;  // order ~1um clearance around all structures

  //  init_default_sector_map();
}

void PHG4CylinderGeom_Spacalv2::ImportParameters(const PHParameters& param)
{
  PHG4CylinderGeom_Spacalv1::ImportParameters(param);

  if (param.exist_int_param("azimuthal_n_sec"))
    azimuthal_n_sec = param.get_int_param("azimuthal_n_sec");
  if (param.exist_double_param("azimuthal_tilt"))
    azimuthal_tilt = param.get_double_param("azimuthal_tilt");
  if (param.exist_int_param("azimuthal_seg_visible"))
    azimuthal_seg_visible = static_cast<bool>(param.get_int_param(
        "azimuthal_seg_visible"));
  if (param.exist_double_param("polar_taper_ratio"))
    polar_taper_ratio = param.get_double_param("polar_taper_ratio");
  if (param.exist_double_param("assembly_spacing"))
    assembly_spacing = param.get_double_param("assembly_spacing");

  return;
}

double
PHG4CylinderGeom_Spacalv2::get_sec_azimuthal_width() const
{
  const double azimuthal_width_base = get_radius() * twopi / (double) (get_azimuthal_n_sec()) - get_assembly_spacing();

  // triggernometry stuff to make a tight connection after tilting

  const double theta1 = get_azimuthal_tilt();
  const double theta2 = pi + (twopi / get_azimuthal_n_sec()) - halfpi - get_azimuthal_tilt();

  return azimuthal_width_base * sin(theta2) / sin(theta1 + theta2);
}

double
PHG4CylinderGeom_Spacalv2::get_half_polar_taper_angle() const
{
  return atan2(get_block_width() * 0.5 * (get_polar_taper_ratio() - 1),
               get_block_depth());
}

int PHG4CylinderGeom_Spacalv2::get_azimuthal_n_sec() const
{
  //  if (config == kNonProjective)
  //    //For kNonProjective geometry, azimuthal_n_sec is calculated, and can not be set externally
  //    return PHG4CylinderGeom_Spacalv1::get_azimuthal_n_sec();
  //  else
  return azimuthal_n_sec;
}

void PHG4CylinderGeom_Spacalv2::set_azimuthal_n_sec(int azimuthalNSec)
{
  if (config == kNonProjective)
  {
    cout
        << "PHG4CylinderGeom_Spacalv2::set_azimuthal_n_sec - Fatal Error - "
           "Spacal is configured as NonProjective geometry. In this case azimuthal_n_sec is calculated, and can not be set externally."
        << endl;
    exit(10);
  }

  azimuthal_n_sec = azimuthalNSec;
  //  init_default_sector_map();
}

bool PHG4CylinderGeom_Spacalv2::is_azimuthal_seg_visible() const
{
  return azimuthal_seg_visible;
}

void PHG4CylinderGeom_Spacalv2::set_azimuthal_seg_visible(bool b)
{
  if (config == kNonProjective)
  {
    cout
        << "PHG4CylinderGeom_Spacalv2::set_azimuthal_seg_visible - Fatal Error - "
           "Spacal is configured as NonProjective geometry. In this case azimuthal_seg_visible is false, and can not be set externally."
        << endl;
    exit(10);
  }

  azimuthal_seg_visible = b;
  ;
}

//! regulated fiber distance in the tapering direction
double
PHG4CylinderGeom_Spacalv2::get_reg_fiber_grid_distance_taper() const
{
  const double mid_plane_width = get_block_width() * ((get_polar_taper_ratio() - 1) * 0.5 + 1) - get_assembly_spacing();

  const int n_grid = floor(mid_plane_width / (get_fiber_distance() * sqrt(3.)));

  return mid_plane_width / n_grid;
}

//! regulated fiber distance in the non-tapering direction
double
PHG4CylinderGeom_Spacalv2::get_reg_fiber_grid_distance_nontaper() const
{
  const double mid_plane_width = get_block_width() - get_assembly_spacing();
  const int n_grid = floor(mid_plane_width / (get_fiber_distance()));
  return mid_plane_width / n_grid;
}
