// $$Id: PHG4CylinderGeom_Spacalv3.cc,v 1.3 2014/08/28 22:18:35 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.3 $$
 * \date $$Date: 2014/08/28 22:18:35 $$
 */

#include "PHG4CylinderGeom_Spacalv3.h"
#include "PHG4CylinderGeom_Spacalv1.h"  // for PHG4CylinderGeom_Spacalv1::...

#include <phparameter/PHParameters.h>

#include <Geant4/G4PhysicalConstants.hh>

#include <CLHEP/Units/SystemOfUnits.h>  // for twopi

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>  // for exit
#include <iostream>
#include <limits>  // std::numeric_limits
#include <map>
#include <sstream>

using namespace std;
using std::make_pair;

PHG4CylinderGeom_Spacalv3::PHG4CylinderGeom_Spacalv3()
{
  SetDefault();
}

PHG4CylinderGeom_Spacalv3::~PHG4CylinderGeom_Spacalv3()
{
  sector_tower_map.erase(sector_tower_map.begin(), sector_tower_map.end());
}

void PHG4CylinderGeom_Spacalv3::identify(std::ostream& os) const
{
  os << "PHG4CylinderGeom_Spacalv3: layer: " << layer
     //
     << ", radius: " << radius
     //
     << ", thickness: " << thickness
     //
     << ", zmin: " << zmin
     //
     << ", zmax: " << zmax << ", num scint: " << nscint << ", num sector: "
     << azimuthal_n_sec << ", unique tower: " << sector_tower_map.size()
     << endl;
}

void PHG4CylinderGeom_Spacalv3::Print(Option_t* opt) const
{
  PHG4CylinderGeom_Spacalv2::Print(opt);

  cout << "\t"
       << "get_sidewall_outer_torr() = " << get_sidewall_outer_torr()
       << endl;
  cout << "\t"
       << "get_sidewall_thickness() = " << get_sidewall_thickness()
       << endl;
  cout << "\t"
       << "get_sidewall_mat() = " << get_sidewall_mat() << endl;
  cout << "\t"
       << "get_max_phi_bin_in_sec() = " << get_max_phi_bin_in_sec()
       << endl;

  subtower_consistency_check();
  cout << "\t"
       << "get_n_subtower_eta() = " << get_n_subtower_eta() << endl;
  cout << "\t"
       << "get_n_subtower_phi() = " << get_n_subtower_phi() << endl;

  cout << "\t"
       << "get_max_lightguide_height() = "
       << get_max_lightguide_height() << endl;
  cout << "\t"
       << "Containing " << sector_tower_map.size()
       << " unique towers per sector." << endl;

  if (get_construction_verbose() >= 2)
    for (const auto& it : sector_tower_map)
    {
      cout << "\t";
      cout << "\t";
      it.second.identify(cout);
    }
}

void PHG4CylinderGeom_Spacalv3::SetDefault()
{
  PHG4CylinderGeom_Spacalv2::SetDefault();

  //  radius = 90.000000;
  //  thickness = 26.130000;
  //  zmin = 149.470000;
  //  zmax = -zmin;
  //  azimuthal_n_sec = 32;
  //  polar_taper_ratio = 1.;
  //  assembly_spacing = 0.002500;
  sidewall_thickness = 0.075000;
  sidewall_outer_torr = 0.030000;
  sidewall_mat = "SS310";
  max_phi_bin_in_sec = 8;
  divider_mat = "G4_AIR";
  divider_width = 14.5;
}

void PHG4CylinderGeom_Spacalv3::ImportParameters(const PHParameters& param)
{
  PHG4CylinderGeom_Spacalv2::ImportParameters(param);

  if (param.exist_double_param("sidewall_thickness"))
    sidewall_thickness = param.get_double_param("sidewall_thickness");
  if (param.exist_double_param("sidewall_outer_torr"))
    sidewall_outer_torr = param.get_double_param("sidewall_outer_torr");
  if (param.exist_string_param("sidewall_mat"))
    sidewall_mat = param.get_string_param("sidewall_mat");
  if (param.exist_int_param("max_phi_bin_in_sec"))
    max_phi_bin_in_sec = param.get_int_param("max_phi_bin_in_sec");
  if (param.exist_string_param("divider_mat"))
    divider_mat = param.get_string_param("divider_mat");
  if (param.exist_double_param("divider_width"))
    divider_width = param.get_double_param("divider_width");

  // load sector_tower_map
  if (param.exist_int_param("sector_tower_map_size"))
  {
    sector_tower_map.clear();

    const int n = param.get_int_param("sector_tower_map_size");

    for (int i = 0; i < n; i++)
    {
      stringstream prefix;
      prefix << "sector_tower_map";
      prefix << "[" << i << "]"
             << ".";

      geom_tower t;
      t.ImportParameters(param, prefix.str());

      sector_tower_map[t.id] = t;
    }
  }

  return;
}

PHG4CylinderGeom_Spacalv3::geom_tower::geom_tower()
  : id(numeric_limits<int>::min())
  ,  //
  pDz(numeric_limits<double>::signaling_NaN())
  ,  //
  pDy1(numeric_limits<double>::signaling_NaN())
  ,  //
  pDx1(numeric_limits<double>::signaling_NaN())
  ,  //
  pDx2(numeric_limits<double>::signaling_NaN())
  ,  //
  pDy2(numeric_limits<double>::signaling_NaN())
  ,  //
  pDx3(numeric_limits<double>::signaling_NaN())
  ,  //
  pDx4(numeric_limits<double>::signaling_NaN())
  ,  //
  pTheta(numeric_limits<double>::signaling_NaN())
  ,  //
  pPhi(numeric_limits<double>::signaling_NaN())
  ,  //
  pAlp1(numeric_limits<double>::signaling_NaN())
  ,  //
  pAlp2(numeric_limits<double>::signaling_NaN())
  ,  //
  pRotationAngleX(numeric_limits<double>::signaling_NaN())
  ,  //
  centralX(numeric_limits<double>::signaling_NaN())
  ,  //
  centralY(numeric_limits<double>::signaling_NaN())
  ,  //
  centralZ(numeric_limits<double>::signaling_NaN())
  ,  //
  ModuleSkinThickness(numeric_limits<double>::signaling_NaN())
  ,  //
  NFiberX(numeric_limits<int>::min())
  ,  //
  NFiberY(numeric_limits<int>::min())
  ,  //
  NSubtowerX(1)
  ,  //
  NSubtowerY(1)
  ,  //
  LightguideHeight(0)
  ,  //
  LightguideTaperRatio(numeric_limits<double>::signaling_NaN())
  ,  //
  LightguideMaterial("PMMA")
{
}

//! fiber layout -> fiber_id
int PHG4CylinderGeom_Spacalv3::geom_tower::compose_fiber_id(int index_x,
                                                            int index_y) const
{
  return NFiberY * index_x + index_y;
}

//! fiber_id -> sub tower ID x: 0 ... NSubtowerX -1
int PHG4CylinderGeom_Spacalv3::geom_tower::get_sub_tower_ID_x(int fiber_id) const
{
  const int index_x = fiber_id / NFiberY;
  assert(index_x < NFiberX);
  assert(index_x >= 0);

  const double sub_tower_width_x = (double) NFiberX / NSubtowerX;
  const int tower_ID_x = (NSubtowerX - 1) - floor(index_x / sub_tower_width_x);  //! x is negative azimuthal direction
  assert(tower_ID_x < NSubtowerX);
  assert(tower_ID_x >= 0);

  return tower_ID_x;
}

//! fiber_id -> sub tower ID y: 0 ... NSubtowerY -1
int PHG4CylinderGeom_Spacalv3::geom_tower::get_sub_tower_ID_y(int fiber_id) const
{
  assert(fiber_id >= 0);
  const int index_y = fiber_id % NFiberY;

  const double sub_tower_width_y = (double) NFiberY / NSubtowerY;

  assert(pRotationAngleX < 0);
  const int tower_ID_y = (NSubtowerY - 1) - floor(index_y / sub_tower_width_y);  //! y is negative polar direction
  assert(tower_ID_y < NSubtowerY);
  assert(tower_ID_y >= 0);

  return tower_ID_y;
}

double
PHG4CylinderGeom_Spacalv3::geom_tower::get_position_fraction_x_in_sub_tower(int fiber_id) const
{
  const int index_x = fiber_id / NFiberY;
  assert(index_x < NFiberX);
  assert(index_x >= 0);

  const double sub_tower_width_x = (double) NFiberX / NSubtowerX;
  const double x_in_sub_tower = (fmod(index_x, sub_tower_width_x) + 0.5) / sub_tower_width_x;  //! x is negative azimuthal direction
  assert(x_in_sub_tower <= 1);
  assert(x_in_sub_tower >= 0);

  return x_in_sub_tower;
}

double
PHG4CylinderGeom_Spacalv3::geom_tower::get_position_fraction_y_in_sub_tower(int fiber_id) const
{
  assert(fiber_id >= 0);
  const int index_y = fiber_id % NFiberY;

  const double sub_tower_width_y = (double) NFiberY / NSubtowerY;

  assert(pRotationAngleX < 0);
  const double y_in_sub_tower = (fmod(index_y, sub_tower_width_y) + 0.5) / sub_tower_width_y;  //! y is negative polar direction
  assert(y_in_sub_tower <= 1);
  assert(y_in_sub_tower >= 0);

  return y_in_sub_tower;
}

void PHG4CylinderGeom_Spacalv3::geom_tower::identify(std::ostream& os) const
{
  os << "PHG4CylinderGeom_Spacalv3::geom_super_tower"
     << "[" << id << "]"
     << " @ <Azimuthal, R, z> = " << centralX << ", " << centralY << ", "
     << centralZ << " cm"
     //
     << " with "
     //
     << "Half length = " << pDz << ", " << pDy1 << ", " << pDx1 << ", " << pDx2
     << ", " << pDy2 << ", " << pDx3 << ", " << pDx4 << ", "
     //
     << "Angles = " << pTheta << ", " << pPhi << ", " << pAlp1 << ", " << pAlp2
     << ", "                              //
     << "Rotation = " << pRotationAngleX  //
     << endl;
}

void PHG4CylinderGeom_Spacalv3::geom_tower::ImportParameters(
    const PHParameters& param, const std::string& param_prefix)
{
  id = param.get_int_param(param_prefix + "id");
  pDz = param.get_double_param(param_prefix + "pDz");

  pDy1 = param.get_double_param(param_prefix + "pDy1");
  pDx1 = param.get_double_param(param_prefix + "pDx1");
  pDx2 = param.get_double_param(param_prefix + "pDx2");
  pDy2 = param.get_double_param(param_prefix + "pDy2");
  pDx3 = param.get_double_param(param_prefix + "pDx3");
  pDx4 = param.get_double_param(param_prefix + "pDx4");

  pTheta = param.get_double_param(param_prefix + "pTheta");
  pPhi = param.get_double_param(param_prefix + "pPhi");
  pAlp1 = param.get_double_param(param_prefix + "pAlp1");
  pAlp2 = param.get_double_param(param_prefix + "pAlp2");

  pRotationAngleX = param.get_double_param(param_prefix + "pRotationAngleX");
  centralX = param.get_double_param(param_prefix + "centralX");
  centralY = param.get_double_param(param_prefix + "centralY");
  centralZ = param.get_double_param(param_prefix + "centralZ");

  ModuleSkinThickness = param.get_double_param(
      param_prefix + "ModuleSkinThickness");
  NFiberX = param.get_int_param(param_prefix + "NFiberX");
  NFiberY = param.get_int_param(param_prefix + "NFiberY");
  NSubtowerX = param.get_int_param(param_prefix + "NSubtowerX");
  NSubtowerY = param.get_int_param(param_prefix + "NSubtowerY");

  LightguideHeight = param.get_double_param(param_prefix + "LightguideHeight");
  LightguideTaperRatio = param.get_double_param(
      param_prefix + "LightguideTaperRatio");
  LightguideMaterial = param.get_string_param(
      param_prefix + "LightguideMaterial");
}

PHG4CylinderGeom_Spacalv3::scint_id_coder::scint_id_coder(int scint_id)
  : scint_ID(scint_id)
{
  sector_ID = (scint_ID >> (kfiber_bit + ktower_bit)) & ((1 << ksector_bit) - 1);
  tower_ID = (scint_ID >> kfiber_bit) & ((1 << ktower_bit) - 1);
  fiber_ID = (scint_ID) & ((1 << kfiber_bit) - 1);
}

PHG4CylinderGeom_Spacalv3::scint_id_coder::scint_id_coder(int sector_id,
                                                          int tower_id, int fiber_id)
  : sector_ID(sector_id)
  , tower_ID(tower_id)
  , fiber_ID(fiber_id)
{
  assert(fiber_ID < (1 << kfiber_bit) and fiber_ID >= 0);
  assert(tower_ID < (1 << ktower_bit) and tower_ID >= 0);
  assert(sector_ID < (1 << ksector_bit) and sector_ID >= 0);

  scint_ID = (((sector_ID << ktower_bit) | tower_ID) << kfiber_bit) | fiber_ID;
}

std::pair<int, int>
PHG4CylinderGeom_Spacalv3::get_tower_z_phi_ID(const int tower_ID,
                                              const int sector_ID) const
{
  // tower_ID to eta/z within a sector
  int z_bin = floor(tower_ID / 10);

  int phi_bin_in_sec = -1;

  if (get_config() == kFullProjective_2DTaper_SameLengthFiberPerTower or get_config() == kFullProjective_2DTaper)
    // colume ID is from -x to +x at the top of the detector, which is reverse of the phi bin direction.
    phi_bin_in_sec = max_phi_bin_in_sec - (tower_ID % 10);
  else if (get_config() == kFullProjective_2DTaper_Tilted_SameLengthFiberPerTower or get_config() == kFullProjective_2DTaper_Tilted)
    phi_bin_in_sec = (tower_ID % 10);

  if (!(phi_bin_in_sec < max_phi_bin_in_sec and phi_bin_in_sec >= 0))
  {
    cout
        << "PHG4CylinderGeom_Spacalv3::get_tower_z_phi_ID - Fatal Error - invalid in put with "
        << "tower_ID = " << tower_ID << ", sector_ID = " << sector_ID << ", phi_bin_in_sec = " << phi_bin_in_sec
        << ". Dump object:" << endl;
    Print();
  }

  assert(phi_bin_in_sec < max_phi_bin_in_sec and phi_bin_in_sec >= 0);

  int phi_bin = sector_ID * max_phi_bin_in_sec + phi_bin_in_sec;

  return make_pair(z_bin, phi_bin);
}

//! get approximate radial position of tower
double PHG4CylinderGeom_Spacalv3::
    get_tower_radial_position(const PHG4CylinderGeom_Spacalv3::geom_tower& tower) const
{
  if (get_config() == kFullProjective_2DTaper_SameLengthFiberPerTower or get_config() == kFullProjective_2DTaper)
    return tower.centralY;
  else if (get_config() == kFullProjective_2DTaper_Tilted_SameLengthFiberPerTower or get_config() == kFullProjective_2DTaper_Tilted)
  {
    const double outter_wall_shift = get_sidewall_thickness() + get_sidewall_outer_torr() + get_assembly_spacing();
    assert(outter_wall_shift >= 0);
    const double tilted_radial_shift = outter_wall_shift / sin(M_PI / get_azimuthal_n_sec());
    assert(tilted_radial_shift >= 0);
    const double tower_radial =  //
        tilted_radial_shift * cos(get_azimuthal_tilt()) +
        get_half_radius() * sin(get_azimuthal_tilt()) * sin(get_azimuthal_tilt()) +
        tower.centralY * cos(get_azimuthal_tilt());

    if (get_construction_verbose() >= 2)
    {
      cout
          << "PHG4CylinderGeom_Spacalv3::get_tower_radial_position - tower radial adjustment: "
             "from "
          << tower.centralY << " to " << tower_radial
          << endl;
    }

    return tower_radial;
  }
  else
  {
    cout
        << "PHG4CylinderGeom_Spacalv3::get_tower_radial_position - ERROR - "
           "unsupported configuration!"
        << endl;
    Print();
    exit(10);
  }
  return NAN;
}

//! check that all towers has consistent sub-tower divider
void PHG4CylinderGeom_Spacalv3::subtower_consistency_check() const
{
  if (sector_tower_map.begin() == sector_tower_map.end()) return;

  for (tower_map_t::const_iterator it = sector_tower_map.begin();
       it != sector_tower_map.end(); ++it)
  {
    assert(get_n_subtower_eta() == it->second.NSubtowerY);
    assert(get_n_subtower_phi() == it->second.NSubtowerX);
  }

  if (get_construction_verbose())
  {
    cout
        << "PHG4CylinderGeom_Spacalv3::subtower_consistency_check - Passed with get_n_subtower_phi() = "
        << get_n_subtower_phi() << " and get_n_subtower_eta()"
        << get_n_subtower_eta() << endl;
  }
}
//! sub-tower divider along the polar direction
int PHG4CylinderGeom_Spacalv3::get_n_subtower_eta() const
{
  if (sector_tower_map.begin() == sector_tower_map.end()) return 0;
  assert(sector_tower_map.begin() != sector_tower_map.end());
  return sector_tower_map.begin()->second.NSubtowerY;
}
//! sub-tower divider along the azimuthal direction
int PHG4CylinderGeom_Spacalv3::get_n_subtower_phi() const
{
  if (sector_tower_map.begin() == sector_tower_map.end()) return 0;
  assert(sector_tower_map.begin() != sector_tower_map.end());
  return sector_tower_map.begin()->second.NSubtowerX;
}

double
PHG4CylinderGeom_Spacalv3::get_max_lightguide_height() const
{
  double max_height = 0;

  for (const auto& it : sector_tower_map)
  {
    const double h = it.second.LightguideHeight;
    max_height = max(max_height, h);
  }

  return max_height;
}

void PHG4CylinderGeom_Spacalv3::load_demo_sector_tower_map1()
{
  cout << "PHG4CylinderGeom_Spacalv3::load_demo_sector_tower_map1 - "
       << "load four example central towers" << endl;

  // Chris Cullen 2D spacal design July 2015
  radius = 90.000000;
  thickness = 26.130000;
  zmin = 149.470000;
  zmax = -zmin;
  azimuthal_n_sec = 32;
  max_phi_bin_in_sec = 8;
  sector_map.clear();
  sector_map[0] = 0;  // only install one sector

  azimuthal_tilt = 0;
  azimuthal_seg_visible = false;
  polar_taper_ratio = 1.;
  assembly_spacing = 0.002500;
  sidewall_thickness = 0.075000;
  sidewall_outer_torr = 0.030000;
  sector_tower_map.clear();

  {
    // tower 1023 based Row/Col = 102/3
    geom_tower geom;
    geom.id = 1023;
    geom.pDz = 6.751948;
    geom.pTheta = 0.038660;
    geom.pPhi = -2.829992;
    geom.pAlp1 = -0.000872;
    geom.pAlp2 = -0.000872;
    geom.pDy1 = 1.121195;
    geom.pDx1 = 1.191868;
    geom.pDx2 = 1.192521;
    geom.pDy2 = 1.281451;
    geom.pDx3 = 1.357679;
    geom.pDx4 = 1.358425;
    geom.centralX = -3.829796;
    geom.centralY = 105.060369;
    geom.centralZ = 3.686651;
    geom.pRotationAngleX = -1.547057;
    geom.ModuleSkinThickness = 0.010000;
    geom.NFiberX = 30;
    geom.NFiberY = 48;
    sector_tower_map[geom.id] = geom;
  }

  {
    // tower 1024 based Row/Col = 102/4
    geom_tower geom;
    geom.id = 1024;
    geom.pDz = 6.751948;
    geom.pTheta = 0.017060;
    geom.pPhi = -2.373142;
    geom.pAlp1 = -0.000290;
    geom.pAlp2 = -0.000290;
    geom.pDy1 = 1.121195;
    geom.pDx1 = 1.190432;
    geom.pDx2 = 1.191082;
    geom.pDy2 = 1.281451;
    geom.pDx3 = 1.356043;
    geom.pDx4 = 1.356787;
    geom.centralX = -1.276086;
    geom.centralY = 105.060369;
    geom.centralZ = 3.686651;
    geom.pRotationAngleX = -1.547057;
    geom.ModuleSkinThickness = 0.010000;
    geom.NFiberX = 30;
    geom.NFiberY = 48;
    sector_tower_map[geom.id] = geom;
  }
}

void PHG4CylinderGeom_Spacalv3::load_demo_sector_tower_map2()
{
  cout << "PHG4CylinderGeom_Spacalv3::load_demo_sector_tower_map2 - "
       << "load one row of example forward towers" << endl;

  // Chris Cullen 2D spacal design July 2015
  radius = 90.000000;
  thickness = 26.130000;
  zmin = 149.470000;
  zmax = -zmin;
  azimuthal_n_sec = 32;
  max_phi_bin_in_sec = 8;
  sector_map.clear();
  sector_map[0] = 0;  // only install one sector
  azimuthal_tilt = 0;
  azimuthal_seg_visible = false;
  polar_taper_ratio = 1.;
  assembly_spacing = 0.002500;
  sidewall_thickness = 0.075000;
  sidewall_outer_torr = 0.030000;
  sector_tower_map.clear();

  {
    // tower 541 based Row/Col = 54/1
    geom_tower geom;
    geom.id = 541;
    geom.pDz = 8.326697;
    geom.pTheta = 0.053578;
    geom.pPhi = -3.015910;
    geom.pAlp1 = 0.068143;
    geom.pAlp2 = 0.068127;
    geom.pDy1 = 1.116997;
    geom.pDx1 = 1.235953;
    geom.pDx2 = 1.214052;
    geom.pDy2 = 1.231781;
    geom.pDx3 = 1.363889;
    geom.pDx4 = 1.339743;
    geom.centralX = -9.019696;
    geom.centralY = 105.782078;
    geom.centralZ = -134.504208;
    geom.pRotationAngleX = -2.482741;
    geom.ModuleSkinThickness = 0.010000;
    geom.NFiberX = 30;
    geom.NFiberY = 48;
    sector_tower_map[geom.id] = geom;
  }
  {
    // tower 542 based Row/Col = 54/2
    geom_tower geom;
    geom.id = 542;
    geom.pDz = 8.326697;
    geom.pTheta = 0.038558;
    geom.pPhi = -2.966436;
    geom.pAlp1 = 0.048643;
    geom.pAlp2 = 0.048632;
    geom.pDy1 = 1.116997;
    geom.pDx1 = 1.234208;
    geom.pDx2 = 1.212398;
    geom.pDy2 = 1.231781;
    geom.pDx3 = 1.361965;
    geom.pDx4 = 1.337918;
    geom.centralX = -6.439665;
    geom.centralY = 105.782078;
    geom.centralZ = -134.504208;
    geom.pRotationAngleX = -2.482741;
    geom.ModuleSkinThickness = 0.010000;
    geom.NFiberX = 30;
    geom.NFiberY = 48;
    sector_tower_map[geom.id] = geom;
  }
  {
    // tower 543 based Row/Col = 54/3
    geom_tower geom;
    geom.id = 543;
    geom.pDz = 8.326697;
    geom.pTheta = 0.023752;
    geom.pPhi = -2.854692;
    geom.pAlp1 = 0.029174;
    geom.pAlp2 = 0.029167;
    geom.pDy1 = 1.116997;
    geom.pDx1 = 1.233047;
    geom.pDx2 = 1.211297;
    geom.pDy2 = 1.231781;
    geom.pDx3 = 1.360684;
    geom.pDx4 = 1.336703;
    geom.centralX = -3.862610;
    geom.centralY = 105.782078;
    geom.centralZ = -134.504208;
    geom.pRotationAngleX = -2.482741;
    geom.ModuleSkinThickness = 0.010000;
    geom.NFiberX = 30;
    geom.NFiberY = 48;
    sector_tower_map[geom.id] = geom;
  }
  {
    // tower 544 based Row/Col = 54/4
    geom_tower geom;
    geom.id = 544;
    geom.pDz = 8.326697;
    geom.pTheta = 0.010142;
    geom.pPhi = -2.416982;
    geom.pAlp1 = 0.009723;
    geom.pAlp2 = 0.009720;
    geom.pDy1 = 1.116997;
    geom.pDx1 = 1.232467;
    geom.pDx2 = 1.210746;
    geom.pDy2 = 1.231781;
    geom.pDx3 = 1.360044;
    geom.pDx4 = 1.336096;
    geom.centralX = -1.287339;
    geom.centralY = 105.782078;
    geom.centralZ = -134.504208;
    geom.pRotationAngleX = -2.482741;
    geom.ModuleSkinThickness = 0.010000;
    geom.NFiberX = 30;
    geom.NFiberY = 48;
    sector_tower_map[geom.id] = geom;
  }
  {
    // tower 545 based Row/Col = 54/5
    geom_tower geom;
    geom.id = 545;
    geom.pDz = 8.326697;
    geom.pTheta = 0.010142;
    geom.pPhi = -0.724610;
    geom.pAlp1 = -0.009723;
    geom.pAlp2 = -0.009720;
    geom.pDy1 = 1.116997;
    geom.pDx1 = 1.232467;
    geom.pDx2 = 1.210746;
    geom.pDy2 = 1.231781;
    geom.pDx3 = 1.360044;
    geom.pDx4 = 1.336096;
    geom.centralX = 1.287339;
    geom.centralY = 105.782078;
    geom.centralZ = -134.504208;
    geom.pRotationAngleX = -2.482741;
    geom.ModuleSkinThickness = 0.010000;
    geom.NFiberX = 30;
    geom.NFiberY = 48;
    sector_tower_map[geom.id] = geom;
  }
  {
    // tower 546 based Row/Col = 54/6
    geom_tower geom;
    geom.id = 546;
    geom.pDz = 8.326697;
    geom.pTheta = 0.023752;
    geom.pPhi = -0.286901;
    geom.pAlp1 = -0.029174;
    geom.pAlp2 = -0.029167;
    geom.pDy1 = 1.116997;
    geom.pDx1 = 1.233047;
    geom.pDx2 = 1.211297;
    geom.pDy2 = 1.231781;
    geom.pDx3 = 1.360684;
    geom.pDx4 = 1.336703;
    geom.centralX = 3.862610;
    geom.centralY = 105.782078;
    geom.centralZ = -134.504208;
    geom.pRotationAngleX = -2.482741;
    geom.ModuleSkinThickness = 0.010000;
    geom.NFiberX = 30;
    geom.NFiberY = 48;
    sector_tower_map[geom.id] = geom;
  }
  {
    // tower 547 based Row/Col = 54/7
    geom_tower geom;
    geom.id = 547;
    geom.pDz = 8.326697;
    geom.pTheta = 0.038558;
    geom.pPhi = -0.175156;
    geom.pAlp1 = -0.048643;
    geom.pAlp2 = -0.048632;
    geom.pDy1 = 1.116997;
    geom.pDx1 = 1.234208;
    geom.pDx2 = 1.212398;
    geom.pDy2 = 1.231781;
    geom.pDx3 = 1.361965;
    geom.pDx4 = 1.337918;
    geom.centralX = 6.439665;
    geom.centralY = 105.782078;
    geom.centralZ = -134.504208;
    geom.pRotationAngleX = -2.482741;
    geom.ModuleSkinThickness = 0.010000;
    geom.NFiberX = 30;
    geom.NFiberY = 48;
    sector_tower_map[geom.id] = geom;
  }
  {
    // tower 548 based Row/Col = 54/8
    geom_tower geom;
    geom.id = 548;
    geom.pDz = 8.326697;
    geom.pTheta = 0.053578;
    geom.pPhi = -0.125683;
    geom.pAlp1 = -0.068143;
    geom.pAlp2 = -0.068127;
    geom.pDy1 = 1.116997;
    geom.pDx1 = 1.235953;
    geom.pDx2 = 1.214052;
    geom.pDy2 = 1.231781;
    geom.pDx3 = 1.363889;
    geom.pDx4 = 1.339743;
    geom.centralX = 9.019696;
    geom.centralY = 105.782078;
    geom.centralZ = -134.504208;
    geom.pRotationAngleX = -2.482741;
    geom.ModuleSkinThickness = 0.010000;
    geom.NFiberX = 30;
    geom.NFiberY = 48;
    sector_tower_map[geom.id] = geom;
  }
}

void PHG4CylinderGeom_Spacalv3::load_demo_sector_tower_map4()
{
  cout << "PHG4CylinderGeom_Spacalv3::load_demo_sector_tower_map4 - "
       << "FermiLab test beam 2014. Need to move to calibration database"
       << endl;

  // From Oleg's documents

  const int ny = 3;
  const int nx = 3;

  const double inch_to_cm = 2.54;

  const double screen_size_y = 2.108 * inch_to_cm;
  const double screen_size_x = 1.049 * inch_to_cm;

  const double z_screen_6_7 = 0.5 * (6.6 + 6.75) * inch_to_cm;
  const double angle_screen_6_7 = 64.78 / 180. * M_PI;
  const double nawrrow_width_x = screen_size_x * sin(angle_screen_6_7);

  const double z_screen_1_2 = 0.5 * (1.35 + 1.6) * inch_to_cm;
  const double angle_screen_1_2 = 90.56 / 180. * M_PI;
  const double wide_width_x = screen_size_x * sin(angle_screen_1_2);

  const double module_length = z_screen_6_7 - z_screen_1_2;
  assert(module_length > 0);

  //tapering, dxwidth/dlength
  const double tapering_ratio = (wide_width_x - nawrrow_width_x) / module_length;
  assert(tapering_ratio < 1);
  assert(tapering_ratio > 0);

  fiber_clading_thickness = 0.003;
  fiber_core_diameter = 0.050 - fiber_clading_thickness * 2;

  assembly_spacing = 0.002500;

  radius = (nawrrow_width_x) / tapering_ratio;
  thickness = module_length * 1.5;  // keep a large torlerence space
  zmin = -(0.5 * ny * screen_size_y + 2 * assembly_spacing * (ny + 1));
  zmax = -zmin;
  azimuthal_n_sec = floor(2 * M_PI / atan(tapering_ratio));
  max_phi_bin_in_sec = 1;

  const double nawrrow_width_x_construction = radius * 2 * tan(M_PI / azimuthal_n_sec) - 2 * assembly_spacing;
  const double wide_width_x_construction = (radius + module_length) * 2 * tan(M_PI / azimuthal_n_sec) - 2 * assembly_spacing;

  cout << "PHG4CylinderGeom_Spacalv3::load_demo_sector_tower_map4 - "

       << "Adjust wide end width by ratio of "
       << wide_width_x_construction / wide_width_x
       << " and narrow end by ratio of "
       << nawrrow_width_x_construction / nawrrow_width_x << endl;

  sector_map.clear();
  for (int sec = 0; sec < nx; ++sec)
  {
    const double rot = twopi / (double) (get_azimuthal_n_sec()) * ((double) (sec) -1);

    sector_map[sec] = rot;
  }

  azimuthal_tilt = 0;
  azimuthal_seg_visible = true;
  polar_taper_ratio = 1.;
  sidewall_thickness = 0;
  sidewall_outer_torr = 0;
  sector_tower_map.clear();

  for (int y = 0; y < ny; y++)
  {
    geom_tower geom;
    geom.id = y;
    geom.pDz = module_length / 2;
    geom.pTheta = 0.;
    geom.pPhi = 0.;
    geom.pAlp1 = 0.;
    geom.pAlp2 = 0.;
    geom.pDy1 = 0.5 * screen_size_y;
    geom.pDx1 = 0.5 * nawrrow_width_x_construction;
    geom.pDx2 = 0.5 * nawrrow_width_x_construction;
    geom.pDy2 = 0.5 * screen_size_y;
    geom.pDx3 = 0.5 * wide_width_x_construction;
    geom.pDx4 = 0.5 * wide_width_x_construction;

    geom.centralX = 0.;
    geom.centralY = module_length * 0.5 + radius + assembly_spacing;
    geom.centralZ = screen_size_y * (y - 1);

    geom.pRotationAngleX = -M_PI / 2.;
    geom.ModuleSkinThickness = 0.010000;
    geom.NFiberX = 30;
    geom.NFiberY = 26 * 2 * 2;
    sector_tower_map[geom.id] = geom;
  }
}
