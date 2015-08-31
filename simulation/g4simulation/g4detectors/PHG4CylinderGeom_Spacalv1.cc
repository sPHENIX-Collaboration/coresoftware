// $$Id: PHG4CylinderGeom_Spacalv1.cc,v 1.3 2014/08/12 03:49:11 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.3 $$
 * \date $$Date: 2014/08/12 03:49:11 $$
 */

#include "PHG4CylinderGeom_Spacalv1.h"

#include <Geant4/globals.hh>
#include <Geant4/G4PhysicalConstants.hh>

#include <cmath>

ClassImp(PHG4CylinderGeom_Spacalv1)

using namespace std;

PHG4CylinderGeom_Spacalv1::PHG4CylinderGeom_Spacalv1()
{
  SetDefault();
}

void
PHG4CylinderGeom_Spacalv1::identify(std::ostream& os) const
{
  os << "PHG4CylinderGeom_Spacalv1: layer: " << layer //
      << ", radius: " << radius //
      << ", thickness: " << thickness //
      << ", zmin: " << zmin //
      << ", zmax: " << zmax << //
      ", num scint: " << nscint

      << endl;
  return;
}

void
PHG4CylinderGeom_Spacalv1::Print(Option_t *) const
{
  identify(cout);

  cout << "Configuration is #" << get_config() << ":" << endl;
  switch (get_config())
    {
  case kNonProjective:
    cout << "fiber always placed radially" << endl;
    break;
  case kProjective_PolarTaper:
    cout
        << "Block constructed with taper in polar direction, non-taper in azimuthal direction. "
        << "The final layout is approximately projective in both azimuthal and polar directions."
        << endl;
    break;
  case kFullProjective_2DTaper:
    cout << "Fully projective spacal with 2D tapered modules" << endl;
    break;
  case kFullProjective_2DTaper_SameLengthFiberPerTower:
    cout
        << "Fully projective spacal with 2D tapered modules. To speed up construction, same-length fiber is used cross one tower"
        << endl;
    break;
  default:
    cout << "PHG4CylinderGeom_Spacalv1::Print - ERROR - unknown configuration #"
        << get_config() << endl;
    break;
    }

  cout << "\t" << "get_max_radius() = " << get_max_radius() << endl;
  cout << "\t" << "get_half_radius() = " << get_half_radius() << endl;
  cout << "\t" << "get_length() = " << get_length() << endl;
  cout << "\t" << "get_*pos() = " << get_xpos() << ", " << get_ypos() << ", "
      << get_zpos() << endl;

  cout << "\t" << "get_azimuthal_n_sec() = " << get_azimuthal_n_sec()
      <<", "<<sector_map.size()<<"/"<< get_azimuthal_n_sec()<<" azimuthal sectors would be filled with SPACAL."<< endl;
  cout << "\t" << "get_azimuthal_distance() = " << get_azimuthal_distance()
      << endl;
  cout << "\t" << "get_z_distance() = " << get_z_distance() << endl;
  cout << "\t" << "get_fiber_outer_r() = " << get_fiber_outer_r() << endl;
  cout << "\t" << "get_fiber_clading_thickness() = "
      << get_fiber_clading_thickness() << endl;
  cout << "\t" << "get_fiber_core_diameter() = " << get_fiber_core_diameter()
      << endl;
  cout << "\t" << "get_fiber_distance() = " << get_fiber_distance() << endl;

  cout << "\t" << "get_absorber_mat() = " << get_absorber_mat() << endl;
  cout << "\t" << "get_fiber_clading_mat() = " << get_fiber_clading_mat()
      << endl;
  cout << "\t" << "get_fiber_core_mat() = " << get_fiber_core_mat() << endl;
  cout << "\t" << "get_calo_step_size() = " << get_calo_step_size() << endl;
  cout << "\t" << "get_fiber_clading_step_size() = "
      << get_fiber_clading_step_size() << endl;
  cout << "\t" << "get_fiber_core_step_size() = " << get_fiber_core_step_size()
      << endl;

  cout << "\t" << "is_virualize_fiber() = " << is_virualize_fiber() << endl;
}

void
PHG4CylinderGeom_Spacalv1::SetDefault()
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

  init_default_sector_map();
}

int
PHG4CylinderGeom_Spacalv1::get_azimuthal_n_sec() const
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
void
PHG4CylinderGeom_Spacalv1::init_default_sector_map()
{
  sector_map.clear();

  for (int sec = 0; sec < get_azimuthal_n_sec(); ++sec)
    {
      const double rot = twopi / (double) (get_azimuthal_n_sec())
          * ((double) (sec));

      sector_map[sec] = rot;
    }
}


