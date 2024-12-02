// $Id: $

/*!
 * \file PHG4CylinderCellGeomSpacalv1.cc
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4CylinderCellGeom_Spacalv1.h"
#include "PHG4CylinderCellDefs.h"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>

void PHG4CylinderCellGeom_Spacalv1::identify(std::ostream& os) const
{
  PHG4CylinderCellGeom::identify(os);

  std::cout << "PHG4CylinderCellGeom_Spacalv1::identify - Tower mapping:" << std::endl;
  for (const tower_z_ID_eta_bin_map_t::value_type& tower_z_ID_eta_bin :
       get_tower_z_ID_eta_bin_map())
  {
    std::cout << "\t"
              << "Tower Z ID[" << tower_z_ID_eta_bin.first
              << "] \t-> Eta Bin " << tower_z_ID_eta_bin.second << std::endl;
  }

  std::cout << "PHG4CylinderCellGeom_Spacalv1::identify - Bin -> z range:" << std::endl;
  for (const bound_map_t::value_type& b : z_bound_map)
  {
    std::cout << "\t"
              << "bin[" << b.first << "] \t-> z = " << b.second.first
              << " - " << b.second.second << std::endl;
  }

  std::cout << "PHG4CylinderCellGeom_Spacalv1::identify - Bin -> eta range:" << std::endl;
  for (const bound_map_t::value_type& b : eta_bound_map)
  {
    std::cout << "\t"
              << "bin[" << b.first << "] \t-> eta = " << b.second.first
              << " - " << b.second.second << std::endl;
  }
  return;
}

void PHG4CylinderCellGeom_Spacalv1::map_consistency_check() const
{
  if ((size_t) nzbins != z_bound_map.size())
  {
    std::cout << "PHG4CylinderCellGeom_Spacalv1::map_consistency_check - "
              << "z_bound_map.size() of " << z_bound_map.size()
              << " in inconsistent with nzbins of " << nzbins << std::endl;
    exit(1);
  }
  if ((size_t) nzbins != eta_bound_map.size())
  {
    std::cout << "PHG4CylinderCellGeom_Spacalv1::map_consistency_check - "
              << "eta_bound_map.size() of " << eta_bound_map.size()
              << " in inconsistent with nzbins of " << nzbins << std::endl;
    exit(1);
  }
  if ((size_t) nzbins < tower_z_ID_eta_bin_map.size())
  {
    std::cout << "PHG4CylinderCellGeom_Spacalv1::map_consistency_check - "
              << "tower_z_ID_eta_bin_map.size() of " << tower_z_ID_eta_bin_map.size()
              << " in inconsistent with nzbins of " << nzbins << std::endl;
    exit(1);
  }
}

void PHG4CylinderCellGeom_Spacalv1::set_zbounds(const int ibin,
                                                const std::pair<double, double>& bounds)
{
  assert(ibin >= 0);
  z_bound_map[ibin] = bounds;
}

void PHG4CylinderCellGeom_Spacalv1::set_etabounds(const int ibin,
                                                  const std::pair<double, double>& bounds)
{
  assert(ibin >= 0);
  eta_bound_map[ibin] = bounds;
}

std::pair<double, double>
PHG4CylinderCellGeom_Spacalv1::get_zbounds(const int ibin) const
{
  map_consistency_check();
  check_binning_method(PHG4CylinderCellDefs::spacalbinning);
  bound_map_t ::const_iterator iter =
      z_bound_map.find(ibin);

  if (iter == z_bound_map.end())
  {
    std::cout
        << "PHG4CylinderCellGeom_Spacalv1::get_zbounds - Fatal Error - Asking for invalid bin in z: "
        << ibin << ". Print of content:" << std::endl;
    identify();
    exit(1);
  }
  return iter->second;
}

std::pair<double, double>
PHG4CylinderCellGeom_Spacalv1::get_etabounds(const int ibin) const
{
  map_consistency_check();
  check_binning_method(PHG4CylinderCellDefs::spacalbinning);

  bound_map_t ::const_iterator iter =
      eta_bound_map.find(ibin);

  if (iter == eta_bound_map.end())
  {
    std::cout
        << "PHG4CylinderCellGeom_Spacalv1::get_etabounds - Fatal Error - Asking for invalid bin in z: "
        << ibin << ". Print of content:" << std::endl;
    identify();
    exit(1);
  }
  return iter->second;
}

int PHG4CylinderCellGeom_Spacalv1::get_zbin(const double /*z*/) const
{
  std::cout << "PHG4CylinderCellGeom_Spacalv1::get_zbin is invalid" << std::endl;
  exit(1);
  return -1;
}

int PHG4CylinderCellGeom_Spacalv1::get_etabin(const double /*eta*/) const
{
  std::cout << "PHG4CylinderCellGeom_Spacalv1::get_etabin is invalid" << std::endl;
  exit(1);
  return -1;
}

double
PHG4CylinderCellGeom_Spacalv1::get_zcenter(const int ibin) const
{
  std::pair<double, double> bound = get_zbounds(ibin);
  return 0.5 * (bound.first + bound.second);
}

double
PHG4CylinderCellGeom_Spacalv1::get_etacenter(const int ibin) const
{
  std::pair<double, double> bound = get_etabounds(ibin);
  return 0.5 * (bound.first + bound.second);
}

int PHG4CylinderCellGeom_Spacalv1::get_etabin_block(const int tower_z_ID) const
{
  map_consistency_check();

  tower_z_ID_eta_bin_map_t::const_iterator iter = tower_z_ID_eta_bin_map.find(tower_z_ID);

  if (iter == tower_z_ID_eta_bin_map.end())
  {
    std::string msg;

    msg = "PHG4CylinderCellGeom_Spacalv1::get_etabin - Fatal Error - can not find tower_z_ID of " + std::to_string(tower_z_ID) + std::string(".");

    throw std::range_error(msg);
  }

  return iter->second;
}
