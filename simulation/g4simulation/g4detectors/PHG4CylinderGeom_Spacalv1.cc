// $$Id: PHG4CylinderGeom_Spacalv1.cc,v 1.3 2014/08/12 03:49:11 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.3 $$
 * \date $$Date: 2014/08/12 03:49:11 $$
 */

#include "PHG4CylinderGeom_Spacalv1.h"

#include <phparameter/PHParameters.h>

#include <Geant4/G4PhysicalConstants.hh>

#include <CLHEP/Units/SystemOfUnits.h>  // for twopi

#include <cmath>
#include <sstream>
#include <utility>  // for pair

PHG4CylinderGeom_Spacalv1::PHG4CylinderGeom_Spacalv1()
{
  SetDefault();
}

void PHG4CylinderGeom_Spacalv1::identify(std::ostream& os) const
{
  os << "PHG4CylinderGeom_Spacalv1: layer: " << layer  //
     << ", radius: " << radius                         //
     << ", thickness: " << thickness                   //
     << ", zmin: " << zmin                             //
     << ", zmax: " << zmax <<                          //
      ", num scint: " << nscint

     << std::endl;
  return;
}

void PHG4CylinderGeom_Spacalv1::Print(Option_t*) const
{
  identify(std::cout);

  std::cout << "Configuration is #" << get_config() << ":" << std::endl;
  switch (get_config())
  {
  case kNonProjective:
    std::cout << "fiber always placed radially" << std::endl;
    break;
  case kFullProjective_2DTaper:
    std::cout << "Fully projective spacal with 2D tapered modules" << std::endl;
    break;
  case kFullProjective_2DTaper_SameLengthFiberPerTower:
    std::cout
        << "Fully projective spacal with 2D tapered modules. To speed up construction, same-length fiber is used cross one tower"
        << std::endl;
    break;
  case kFullProjective_2DTaper_Tilted:
    std::cout << "Fully projective spacal with 2D tapered modules and  allow azimuthal tilts" << std::endl;
    break;
  case kFullProjective_2DTaper_Tilted_SameLengthFiberPerTower:
    std::cout
        << "Fully projective spacal with 2D tapered modules and  allow azimuthal tilts. To speed up construction, same-length fiber is used cross one tower"
        << std::endl;
    break;
  default:
    std::cout << "PHG4CylinderGeom_Spacalv1::Print - ERROR - unknown configuration #"
              << get_config() << std::endl;
    break;
  }

  std::cout << "\t"
            << "get_max_radius() = " << get_max_radius() << std::endl;
  std::cout << "\t"
            << "get_half_radius() = " << get_half_radius() << std::endl;
  std::cout << "\t"
            << "get_length() = " << get_length() << std::endl;
  std::cout << "\t"
            << "get_*pos() = " << get_xpos() << ", " << get_ypos() << ", "
            << get_zpos() << std::endl;

  std::cout << "\t"
            << "get_azimuthal_n_sec() = " << get_azimuthal_n_sec() << ", "
            << sector_map.size() << "/" << get_azimuthal_n_sec()
            << " azimuthal sectors would be filled with SPACAL." << std::endl;
  std::cout << "\t"
            << "get_azimuthal_distance() = " << get_azimuthal_distance()
            << std::endl;
  std::cout << "\t"
            << "get_z_distance() = " << get_z_distance() << std::endl;
  std::cout << "\t"
            << "get_fiber_outer_r() = " << get_fiber_outer_r() << std::endl;
  std::cout << "\t"
            << "get_fiber_clading_thickness() = "
            << get_fiber_clading_thickness() << std::endl;
  std::cout << "\t"
            << "get_fiber_core_diameter() = " << get_fiber_core_diameter()
            << std::endl;
  std::cout << "\t"
            << "get_fiber_distance() = " << get_fiber_distance() << std::endl;

  std::cout << "\t"
            << "get_absorber_mat() = " << get_absorber_mat() << std::endl;
  std::cout << "\t"
            << "get_fiber_clading_mat() = " << get_fiber_clading_mat()
            << std::endl;
  std::cout << "\t"
            << "get_fiber_core_mat() = " << get_fiber_core_mat() << std::endl;
  //  std::cout << "\t" << "get_calo_step_size() = " << get_calo_step_size() << std::endl;
  //  std::cout << "\t" << "get_fiber_clading_step_size() = "
  //      << get_fiber_clading_step_size() << std::endl;
  std::cout << "\t"
            << "get_fiber_core_step_size() = " << get_fiber_core_step_size()
            << std::endl;

  std::cout << "\t"
            << "is_virualize_fiber() = " << is_virualize_fiber() << std::endl;
  std::cout << "\t"
            << "get_construction_verbose() = " << get_construction_verbose()
            << std::endl;

  if (get_construction_verbose() >= 2)
  {
    std::cout << "\t"
              << "Containing " << sector_map.size()
              << " sector with rotation specified:" << std::endl;
    for (auto it : sector_map)
    {
      std::cout << "\t"
                << "\t"
                << "sector_map[" << it.first << "] = " << it.second
                << std::endl;
    }
  }
}

void PHG4CylinderGeom_Spacalv1::SetDefault()
{
  config = kNonProjective;

  layer = 0;
  radius = 95;
  thickness = 16.6;
  zmin = -143;
  zmax = -zmin;
  nscint = 0;

  absorber_mat = "Spacal_W_Epoxy";
  fiber_clading_mat = "PMMA";
  fiber_core_mat = "G4_POLYSTYRENE";

  xpos = 0;
  ypos = 0;
  zpos = 0;

  fiber_clading_thickness = 0.003 / 2;
  fiber_core_diameter = 0.047 - fiber_clading_thickness * 2;
  fiber_distance = 0.1;

  virualize_fiber = false;
  construction_verbose = 0;

  //  init_default_sector_map();
}

void PHG4CylinderGeom_Spacalv1::ImportParameters(const PHParameters& param)
{
  PHG4CylinderGeomv2::ImportParameters(param);

  if (param.exist_string_param("absorber_mat"))
  {
    absorber_mat = param.get_string_param("absorber_mat");
  }
  if (param.exist_string_param("fiber_core_mat"))
  {
    fiber_core_mat = param.get_string_param("fiber_core_mat");
  }
  if (param.exist_string_param("fiber_clading_mat"))
  {
    fiber_clading_mat = param.get_string_param("fiber_clading_mat");
  }
  if (param.exist_double_param("xpos"))
  {
    xpos = param.get_double_param("xpos");
  }
  if (param.exist_double_param("ypos"))
  {
    ypos = param.get_double_param("ypos");
  }
  if (param.exist_double_param("zpos"))
  {
    zpos = param.get_double_param("zpos");
  }
  if (param.exist_double_param("fiber_core_diameter"))
  {
    fiber_core_diameter = param.get_double_param("fiber_core_diameter");
  }
  if (param.exist_double_param("fiber_clading_thickness"))
  {
    fiber_clading_thickness = param.get_double_param("fiber_clading_thickness");
  }
  if (param.exist_double_param("fiber_distance"))
  {
    fiber_distance = param.get_double_param("fiber_distance");
  }
  if (param.exist_int_param("config"))
  {
    config = static_cast<config_t>(param.get_int_param("config"));
  }
  if (param.exist_int_param("virualize_fiber"))
  {
    virualize_fiber = static_cast<bool>(param.get_int_param("virualize_fiber"));
  }
  if (param.exist_int_param("construction_verbose"))
  {
    construction_verbose = param.get_int_param("construction_verbose");
  }

  // init_default_sector_map if instructed to do so
  if (param.exist_int_param("init_default_sector_map"))
  {
    if (param.get_int_param("init_default_sector_map"))
    {
      init_default_sector_map();
    }
  }

  // load sector_map if specified. Over write init_default_sector_map if both presents
  if (param.exist_int_param("sector_map_size"))
  {
    sector_map.clear();

    const int n = param.get_int_param("sector_map_size");

    for (int i = 0; i < n; i++)
    {
      std::stringstream prefix;
      prefix << "sector_map";
      prefix << "[" << i << "]"
             << ".";

      const int id = param.get_int_param(prefix.str() + "id");
      const double rotation = param.get_double_param(
          prefix.str() + "rotation");

      sector_map[id] = rotation;
    }
  }

  return;
}

int PHG4CylinderGeom_Spacalv1::get_azimuthal_n_sec() const
{
  return std::floor(
      get_half_radius() * twopi / (get_fiber_distance() * sqrt(3.)));
}

double
PHG4CylinderGeom_Spacalv1::get_azimuthal_distance() const
{
  return get_half_radius() * twopi / (double) (get_azimuthal_n_sec());
}

double
PHG4CylinderGeom_Spacalv1::get_z_distance() const
{
  return get_fiber_distance() / 2.;
}

//! load a default map that populate all the sectors
void PHG4CylinderGeom_Spacalv1::init_default_sector_map()
{
  sector_map.clear();

  for (int sec = 0; sec < get_azimuthal_n_sec(); ++sec)
  {
    const double rot = twopi / (double) (get_azimuthal_n_sec()) * ((double) (sec));

    sector_map[sec] = rot;
  }
}
