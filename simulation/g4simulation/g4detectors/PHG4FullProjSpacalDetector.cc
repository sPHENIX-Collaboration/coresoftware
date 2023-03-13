
/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.3 $$
 * \date $$Date: 2015/02/10 15:39:26 $$
 */
#include "PHG4FullProjSpacalDetector.h"

#include "PHG4SpacalDisplayAction.h"

#include <g4gdml/PHG4GDMLConfig.hh>

#include <phool/recoConsts.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Exception.hh>          // for G4Exception, G4ExceptionD
#include <Geant4/G4ExceptionSeverity.hh>  // for FatalException
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Transform3D.hh>
#include <Geant4/G4Trap.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>  // for G4double
#include <Geant4/G4Vector3D.hh>

#include <TSystem.h>

#include <boost/foreach.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>  // for operator<<, basic_ostream
#include <map>       // for allocator, map<>::value_type
#include <numeric>   // std::accumulate
#include <sstream>
#include <string>  // std::string, std::to_string
#include <vector>  // for vector

class G4Material;
class PHCompositeNode;

//_______________________________________________________________
//note this inactive thickness is ~1.5% of a radiation length
PHG4FullProjSpacalDetector::PHG4FullProjSpacalDetector(PHG4Subsystem* subsys, PHCompositeNode* Node,
                                                       const std::string& dnam, PHParameters* parameters, const int lyr)
  : PHG4SpacalDetector(subsys, Node, dnam, parameters, lyr, false)
{
  assert(_geom == nullptr);

  _geom = new SpacalGeom_t();
  if (_geom == nullptr)
  {
    std::cout
        << "PHG4FullProjSpacalDetector::Constructor - Fatal Error - invalid geometry object!"
        << std::endl;
    gSystem->Exit(1);
  }

  //this class loads Chris Cullen 2D spacal design July 2015 by default.
  // this step is deprecated now
  // get_geom_v3()->load_demo_sector_tower_map_2015_Chris_Cullen_2D_spacal();

  assert(parameters);
  get_geom_v3()->ImportParameters(*parameters);

  //  std::cout <<"PHG4FullProjSpacalDetector::Constructor -  get_geom_v3()->Print();"<<std::endl;
  //  get_geom_v3()->Print();
}

//_______________________________________________________________
void PHG4FullProjSpacalDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if (get_geom_v3()->get_construction_verbose() >= 1)
  {
    std::cout << "PHG4FullProjSpacalDetector::Construct::" << GetName()
              << " - start with PHG4SpacalDetector::Construct()." << std::endl;
  }

  PHG4SpacalDetector::ConstructMe(logicWorld);

  if (get_geom_v3()->get_construction_verbose() >= 1)
  {
    std::cout << "PHG4FullProjSpacalDetector::Construct::" << GetName()
              << " - Completed." << std::endl;
  }
}

std::pair<G4LogicalVolume*, G4Transform3D>
PHG4FullProjSpacalDetector::Construct_AzimuthalSeg()
{
  if (!(get_geom_v3()->get_azimuthal_n_sec() > 4))
  {
    std::cout << "azimuthal_n_sec <= 4: " << get_geom_v3()->get_azimuthal_n_sec() << std::endl;
    gSystem->Exit(1);
  }

  G4Tubs* sec_solid = new G4Tubs(G4String(GetName() + std::string("_sec")),
                                 get_geom_v3()->get_radius() * cm, get_geom_v3()->get_max_radius() * cm,
                                 get_geom_v3()->get_length() * cm / 2.0,
                                 halfpi - pi / get_geom_v3()->get_azimuthal_n_sec(),
                                 twopi / get_geom_v3()->get_azimuthal_n_sec());

  recoConsts* rc = recoConsts::instance();
  G4Material* cylinder_mat = GetDetectorMaterial(rc->get_StringFlag("WorldMaterial"));
  assert(cylinder_mat);

  G4LogicalVolume* sec_logic = new G4LogicalVolume(sec_solid, cylinder_mat,
                                                   G4String(G4String(GetName() + std::string("_sec"))), nullptr, nullptr);

  GetDisplayAction()->AddVolume(sec_logic, "Sector");

  // construct walls

  G4Material* wall_mat = GetDetectorMaterial(get_geom_v3()->get_sidewall_mat());
  assert(wall_mat);

  if (get_geom_v3()->get_sidewall_thickness() > 0)
  {
    // end walls
    if (get_geom_v3()->get_construction_verbose() >= 1)
    {
      std::cout << "PHG4FullProjSpacalDetector::Construct_AzimuthalSeg::" << GetName()
                << " - construct end walls." << std::endl;
    }
    G4Tubs* wall_solid = new G4Tubs(G4String(GetName() + std::string("_EndWall")),
                                    get_geom_v3()->get_radius() * cm + get_geom_v3()->get_sidewall_outer_torr() * cm,
                                    get_geom_v3()->get_max_radius() * cm - get_geom_v3()->get_sidewall_outer_torr() * cm,
                                    get_geom_v3()->get_sidewall_thickness() * cm / 2.0,
                                    halfpi - pi / get_geom_v3()->get_azimuthal_n_sec(),
                                    twopi / get_geom_v3()->get_azimuthal_n_sec());

    G4LogicalVolume* wall_logic = new G4LogicalVolume(wall_solid, wall_mat,
                                                      G4String(G4String(GetName() + std::string("_EndWall"))), nullptr, nullptr,
                                                      nullptr);
    GetDisplayAction()->AddVolume(wall_logic, "WallProj");

    using z_locations_t = std::map<int, double>;
    z_locations_t z_locations;
    z_locations[1000] = get_geom_v3()->get_sidewall_thickness() * cm / 2.0 + get_geom_v3()->get_assembly_spacing() * cm;
    z_locations[1001] = get_geom_v3()->get_length() * cm / 2.0 - (get_geom_v3()->get_sidewall_thickness() * cm / 2.0 + get_geom_v3()->get_assembly_spacing() * cm);
    z_locations[1100] = -(get_geom_v3()->get_sidewall_thickness() * cm / 2.0 + get_geom_v3()->get_assembly_spacing() * cm);
    z_locations[1101] = -(get_geom_v3()->get_length() * cm / 2.0 - (get_geom_v3()->get_sidewall_thickness() * cm / 2.0 + get_geom_v3()->get_assembly_spacing() * cm));

    BOOST_FOREACH (z_locations_t::value_type& val, z_locations)
    {
      if (get_geom_v3()->get_construction_verbose() >= 2)
        std::cout << "PHG4FullProjSpacalDetector::Construct_AzimuthalSeg::"
                  << GetName() << " - constructed End Wall ID " << val.first
                  << " @ Z = " << val.second << std::endl;

      G4Transform3D wall_trans = G4TranslateZ3D(val.second);

      G4PVPlacement* wall_phys = new G4PVPlacement(wall_trans, wall_logic,
                                                   G4String(GetName()) + G4String("_EndWall"), sec_logic,
                                                   false, val.first, OverlapCheck());

      calo_vol[wall_phys] = val.first;
      assert(gdml_config);
      gdml_config->exclude_physical_vol(wall_phys);
    }
  }

  if (get_geom_v3()->get_sidewall_thickness() > 0)
  {
    // side walls
    if (get_geom_v3()->get_construction_verbose() >= 1)
    {
      std::cout << "PHG4FullProjSpacalDetector::Construct_AzimuthalSeg::" << GetName()
                << " - construct side walls." << std::endl;
    }
    G4Box* wall_solid = new G4Box(G4String(GetName() + std::string("_SideWall")),
                                  get_geom_v3()->get_sidewall_thickness() * cm / 2.0,
                                  get_geom_v3()->get_thickness() * cm / 2. - 2 * get_geom_v3()->get_sidewall_outer_torr() * cm,
                                  (get_geom_v3()->get_length() / 2. - 2 * (get_geom_v3()->get_sidewall_thickness() + 2. * get_geom_v3()->get_assembly_spacing())) * cm * .5);

    G4LogicalVolume* wall_logic = new G4LogicalVolume(wall_solid, wall_mat,
                                                      G4String(G4String(GetName() + std::string("_SideWall"))), nullptr, nullptr,
                                                      nullptr);
    GetDisplayAction()->AddVolume(wall_logic, "WallProj");

    using sign_t = std::map<int, std::pair<int, int>>;
    sign_t signs;
    signs[2000] = std::make_pair(+1, +1);
    signs[2001] = std::make_pair(+1, -1);
    signs[2100] = std::make_pair(-1, +1);
    signs[2101] = std::make_pair(-1, -1);

    BOOST_FOREACH (sign_t::value_type& val, signs)
    {
      const int sign_z = val.second.first;
      const int sign_azimuth = val.second.second;

      if (get_geom_v3()->get_construction_verbose() >= 2)
        std::cout << "PHG4FullProjSpacalDetector::Construct_AzimuthalSeg::"
                  << GetName() << " - constructed Side Wall ID " << val.first
                  << " with"
                  << " Shift X = "
                  << sign_azimuth * (get_geom_v3()->get_sidewall_thickness() * cm / 2.0 + get_geom_v3()->get_sidewall_outer_torr() * cm)
                  << " Rotation Z = "
                  << sign_azimuth * pi / get_geom_v3()->get_azimuthal_n_sec()
                  << " Shift Z = " << sign_z * (get_geom_v3()->get_length() * cm / 4)
                  << std::endl;

      G4Transform3D wall_trans = G4RotateZ3D(
                                     sign_azimuth * pi / get_geom_v3()->get_azimuthal_n_sec()) *
                                 G4TranslateZ3D(sign_z * (get_geom_v3()->get_length() * cm / 4)) * G4TranslateY3D(get_geom_v3()->get_radius() * cm + get_geom_v3()->get_thickness() * cm / 2.) * G4TranslateX3D(sign_azimuth * (get_geom_v3()->get_sidewall_thickness() * cm / 2.0 + get_geom_v3()->get_sidewall_outer_torr() * cm));

      G4PVPlacement* wall_phys = new G4PVPlacement(wall_trans, wall_logic,
                                                   G4String(GetName()) + G4String("_EndWall"), sec_logic,
                                                   false, val.first, OverlapCheck());

      calo_vol[wall_phys] = val.first;

      assert(gdml_config);
      gdml_config->exclude_physical_vol(wall_phys);
    }
  }

  // construct towers

  BOOST_FOREACH (const SpacalGeom_t::tower_map_t::value_type& val, get_geom_v3()->get_sector_tower_map())
  {
    const SpacalGeom_t::geom_tower& g_tower = val.second;
    G4LogicalVolume* LV_tower = Construct_Tower(g_tower);

    G4Transform3D block_trans = G4TranslateX3D(g_tower.centralX * cm) * G4TranslateY3D(g_tower.centralY * cm) * G4TranslateZ3D(g_tower.centralZ * cm) * G4RotateX3D(g_tower.pRotationAngleX * rad);

    const bool overlapcheck_block = OverlapCheck() and (get_geom_v3()->get_construction_verbose() >= 2);

    G4PVPlacement* block_phys = new G4PVPlacement(block_trans, LV_tower,
                                                  G4String(GetName()) + G4String("_Tower"), sec_logic, false,
                                                  g_tower.id, overlapcheck_block);
    block_vol[block_phys] = g_tower.id;

    assert(gdml_config);
    gdml_config->exclude_physical_vol(block_phys);
  }

  std::cout << "PHG4FullProjSpacalDetector::Construct_AzimuthalSeg::" << GetName()
            << " - constructed " << get_geom_v3()->get_sector_tower_map().size()
            << " unique towers" << std::endl;

  return std::make_pair(sec_logic, G4Transform3D::Identity);
}

//! Fully projective spacal with 2D tapered modules. To speed up construction, same-length fiber is used cross one tower
int PHG4FullProjSpacalDetector::Construct_Fibers_SameLengthFiberPerTower(
    const PHG4FullProjSpacalDetector::SpacalGeom_t::geom_tower& g_tower,
    G4LogicalVolume* LV_tower)
{
  // construct fibers

  // first check out the fibers geometry

  using fiber_par_map = std::map<int, std::pair<G4Vector3D, G4Vector3D>>;
  fiber_par_map fiber_par;
  G4double min_fiber_length = g_tower.pDz * cm * 4;

  G4Vector3D v_zshift = G4Vector3D(tan(g_tower.pTheta) * cos(g_tower.pPhi),
                                   tan(g_tower.pTheta) * sin(g_tower.pPhi), 1) *
                        g_tower.pDz;
  //  int fiber_ID = 0;
  for (int ix = 0; ix < g_tower.NFiberX; ix++)
  //  int ix = 0;
  {
    const double weighted_ix = static_cast<double>(ix) / (g_tower.NFiberX - 1.);

    const double weighted_pDx1 = (g_tower.pDx1 - g_tower.ModuleSkinThickness - get_geom_v3()->get_fiber_outer_r()) * (weighted_ix * 2 - 1);
    const double weighted_pDx2 = (g_tower.pDx2 - g_tower.ModuleSkinThickness - get_geom_v3()->get_fiber_outer_r()) * (weighted_ix * 2 - 1);

    const double weighted_pDx3 = (g_tower.pDx3 - g_tower.ModuleSkinThickness - get_geom_v3()->get_fiber_outer_r()) * (weighted_ix * 2 - 1);
    const double weighted_pDx4 = (g_tower.pDx4 - g_tower.ModuleSkinThickness - get_geom_v3()->get_fiber_outer_r()) * (weighted_ix * 2 - 1);

    for (int iy = 0; iy < g_tower.NFiberY; iy++)
    //        int iy = 0;
    {
      if ((ix + iy) % 2 == 1)
        continue;  // make a triangle pattern

      const double weighted_iy = static_cast<double>(iy) / (g_tower.NFiberY - 1.);

      const double weighted_pDy1 = (g_tower.pDy1 - g_tower.ModuleSkinThickness - get_geom_v3()->get_fiber_outer_r()) * (weighted_iy * 2 - 1);
      const double weighted_pDy2 = (g_tower.pDy2 - g_tower.ModuleSkinThickness - get_geom_v3()->get_fiber_outer_r()) * (weighted_iy * 2 - 1);

      const double weighted_pDx12 = weighted_pDx1 * (1 - weighted_iy) + weighted_pDx2 * (weighted_iy) + weighted_pDy1 * tan(g_tower.pAlp1);
      const double weighted_pDx34 = weighted_pDx3 * (1 - weighted_iy) + weighted_pDx4 * (weighted_iy) + weighted_pDy1 * tan(g_tower.pAlp2);

      G4Vector3D v1 = G4Vector3D(weighted_pDx12, weighted_pDy1, 0) - v_zshift;
      G4Vector3D v2 = G4Vector3D(weighted_pDx34, weighted_pDy2, 0) + v_zshift;

      G4Vector3D vector_fiber = (v2 - v1);
      vector_fiber *= (vector_fiber.mag() - get_geom_v3()->get_fiber_outer_r()) / vector_fiber.mag();  // shrink by fiber boundary protection
      G4Vector3D center_fiber = (v2 + v1) / 2;

      // convert to Geant4 units
      vector_fiber *= cm;
      center_fiber *= cm;

      const int fiber_ID = g_tower.compose_fiber_id(ix, iy);
      fiber_par[fiber_ID] = std::make_pair(vector_fiber,
                                           center_fiber);

      const G4double fiber_length = vector_fiber.mag();

      min_fiber_length = std::min(fiber_length, min_fiber_length);

      //          ++fiber_ID;
    }
  }

  int fiber_count = 0;

  const G4double fiber_length = min_fiber_length;
  std::vector<G4double> fiber_cut;

  std::stringstream ss;
  ss << std::string("_Tower") << g_tower.id;
  G4LogicalVolume* fiber_logic = Construct_Fiber(fiber_length, ss.str());

  BOOST_FOREACH (const fiber_par_map::value_type& val, fiber_par)
  {
    const int fiber_ID = val.first;
    G4Vector3D vector_fiber = val.second.first;
    G4Vector3D center_fiber = val.second.second;
    const G4double optimal_fiber_length = vector_fiber.mag();

    const G4Vector3D v1 = center_fiber - 0.5 * vector_fiber;

    // keep a statistics
    assert(optimal_fiber_length - fiber_length >= 0);
    fiber_cut.push_back(optimal_fiber_length - fiber_length);

    center_fiber += (fiber_length / optimal_fiber_length - 1) * 0.5 * vector_fiber;
    vector_fiber *= fiber_length / optimal_fiber_length;

    //      const G4Vector3D v1_new = center_fiber - 0.5 *vector_fiber;

    if (get_geom_v3()->get_construction_verbose() >= 3)
      std::cout << "PHG4FullProjSpacalDetector::Construct_Fibers_SameLengthFiberPerTower::" << GetName()
                << " - constructed fiber " << fiber_ID << ss.str()  //
                << ", Length = " << optimal_fiber_length << "-"
                << (optimal_fiber_length - fiber_length) << "mm, "  //
                << "x = " << center_fiber.x() << "mm, "             //
                << "y = " << center_fiber.y() << "mm, "             //
                << "z = " << center_fiber.z() << "mm, "             //
                << "vx = " << vector_fiber.x() << "mm, "            //
                << "vy = " << vector_fiber.y() << "mm, "            //
                << "vz = " << vector_fiber.z() << "mm, "            //
                << std::endl;

    const G4double rotation_angle = G4Vector3D(0, 0, 1).angle(vector_fiber);
    const G4Vector3D rotation_axis =
        rotation_angle == 0 ? G4Vector3D(1, 0, 0) : G4Vector3D(0, 0, 1).cross(vector_fiber);

    G4Transform3D fiber_place(
        G4Translate3D(center_fiber.x(), center_fiber.y(), center_fiber.z()) * G4Rotate3D(rotation_angle, rotation_axis));

    std::stringstream name;
    name << GetName() + std::string("_Tower") << g_tower.id << "_fiber"
         << ss.str();

    const bool overlapcheck_fiber = OverlapCheck() and (get_geom_v3()->get_construction_verbose() >= 3);
    G4PVPlacement* fiber_physi = new G4PVPlacement(fiber_place, fiber_logic,
                                                   G4String(name.str()), LV_tower, false, fiber_ID,
                                                   overlapcheck_fiber);
    fiber_vol[fiber_physi] = fiber_ID;

    assert(gdml_config);
    gdml_config->exclude_physical_vol(fiber_physi);

    fiber_count++;
  }

  if (get_geom_v3()->get_construction_verbose() >= 2)
    std::cout
        << "PHG4FullProjSpacalDetector::Construct_Fibers_SameLengthFiberPerTower::"
        << GetName() << " - constructed tower ID " << g_tower.id << " with "
        << fiber_count << " fibers. Average fiber length cut = "
        << accumulate(fiber_cut.begin(), fiber_cut.end(), 0.0) / fiber_cut.size() << " mm" << std::endl;

  return fiber_count;
}

//! a block along z axis built with G4Trd that is slightly tapered in x dimension
int PHG4FullProjSpacalDetector::Construct_Fibers(
    const PHG4FullProjSpacalDetector::SpacalGeom_t::geom_tower& g_tower,
    G4LogicalVolume* LV_tower)
{
  G4Vector3D v_zshift = G4Vector3D(tan(g_tower.pTheta) * cos(g_tower.pPhi),
                                   tan(g_tower.pTheta) * sin(g_tower.pPhi), 1) *
                        g_tower.pDz;
  int fiber_cnt = 0;
  for (int ix = 0; ix < g_tower.NFiberX; ix++)
  {
    const double weighted_ix = static_cast<double>(ix) / (g_tower.NFiberX - 1.);

    const double weighted_pDx1 = (g_tower.pDx1 - g_tower.ModuleSkinThickness - get_geom_v3()->get_fiber_outer_r()) * (weighted_ix * 2 - 1);
    const double weighted_pDx2 = (g_tower.pDx2 - g_tower.ModuleSkinThickness - get_geom_v3()->get_fiber_outer_r()) * (weighted_ix * 2 - 1);

    const double weighted_pDx3 = (g_tower.pDx3 - g_tower.ModuleSkinThickness - get_geom_v3()->get_fiber_outer_r()) * (weighted_ix * 2 - 1);
    const double weighted_pDx4 = (g_tower.pDx4 - g_tower.ModuleSkinThickness - get_geom_v3()->get_fiber_outer_r()) * (weighted_ix * 2 - 1);

    for (int iy = 0; iy < g_tower.NFiberY; iy++)
    {
      if ((ix + iy) % 2 == 1)
        continue;  // make a triangle pattern
      const int fiber_ID = g_tower.compose_fiber_id(ix, iy);

      const double weighted_iy = static_cast<double>(iy) / (g_tower.NFiberY - 1.);

      const double weighted_pDy1 = (g_tower.pDy1 - g_tower.ModuleSkinThickness - get_geom_v3()->get_fiber_outer_r()) * (weighted_iy * 2 - 1);
      const double weighted_pDy2 = (g_tower.pDy2 - g_tower.ModuleSkinThickness - get_geom_v3()->get_fiber_outer_r()) * (weighted_iy * 2 - 1);

      const double weighted_pDx12 = weighted_pDx1 * (1 - weighted_iy) + weighted_pDx2 * (weighted_iy) + weighted_pDy1 * tan(g_tower.pAlp1);
      const double weighted_pDx34 = weighted_pDx3 * (1 - weighted_iy) + weighted_pDx4 * (weighted_iy) + weighted_pDy1 * tan(g_tower.pAlp2);

      G4Vector3D v1 = G4Vector3D(weighted_pDx12, weighted_pDy1, 0) - v_zshift;
      G4Vector3D v2 = G4Vector3D(weighted_pDx34, weighted_pDy2, 0) + v_zshift;

      G4Vector3D vector_fiber = (v2 - v1);
      vector_fiber *= (vector_fiber.mag() - get_geom_v3()->get_fiber_outer_r()) / vector_fiber.mag();  // shrink by fiber boundary protection
      G4Vector3D center_fiber = (v2 + v1) / 2;

      // convert to Geant4 units
      vector_fiber *= cm;
      center_fiber *= cm;

      const G4double fiber_length = vector_fiber.mag();

      std::stringstream ss;
      ss << std::string("_Tower") << g_tower.id;
      ss << "_x" << ix;
      ss << "_y" << iy;
      G4LogicalVolume* fiber_logic = Construct_Fiber(fiber_length,
                                                     ss.str());

      if (get_geom_v3()->get_construction_verbose() >= 3)
        std::cout << "PHG4FullProjSpacalDetector::Construct_Fibers::" << GetName()
                  << " - constructed fiber " << fiber_ID << ss.str()  //
                  << ", Length = " << fiber_length << "mm, "          //
                  << "x = " << center_fiber.x() << "mm, "             //
                  << "y = " << center_fiber.y() << "mm, "             //
                  << "z = " << center_fiber.z() << "mm, "             //
                  << "vx = " << vector_fiber.x() << "mm, "            //
                  << "vy = " << vector_fiber.y() << "mm, "            //
                  << "vz = " << vector_fiber.z() << "mm, "            //
                  << std::endl;

      const G4double rotation_angle = G4Vector3D(0, 0, 1).angle(
          vector_fiber);
      const G4Vector3D rotation_axis =
          rotation_angle == 0 ? G4Vector3D(1, 0, 0) : G4Vector3D(0, 0, 1).cross(vector_fiber);

      G4Transform3D fiber_place(
          G4Translate3D(center_fiber.x(), center_fiber.y(),
                        center_fiber.z()) *
          G4Rotate3D(rotation_angle, rotation_axis));

      std::stringstream name;
      name << GetName() + std::string("_Tower") << g_tower.id << "_fiber"
           << ss.str();

      const bool overlapcheck_fiber = OverlapCheck() and (get_geom_v3()->get_construction_verbose() >= 3);
      G4PVPlacement* fiber_physi = new G4PVPlacement(fiber_place,
                                                     fiber_logic, G4String(name.str()), LV_tower, false,
                                                     fiber_ID, overlapcheck_fiber);
      fiber_vol[fiber_physi] = fiber_ID;

      assert(gdml_config);
      gdml_config->exclude_physical_vol(fiber_physi);

      ++fiber_cnt;
    }
  }

  if (get_geom_v3()->get_construction_verbose() >= 3)
    std::cout << "PHG4FullProjSpacalDetector::Construct_Fibers::" << GetName()
              << " - constructed tower ID " << g_tower.id << " with " << fiber_cnt
              << " fibers" << std::endl;

  return fiber_cnt;
}

//! a block along z axis built with G4Trd that is slightly tapered in x dimension
G4LogicalVolume*
PHG4FullProjSpacalDetector::Construct_Tower(
    const PHG4FullProjSpacalDetector::SpacalGeom_t::geom_tower& g_tower)
{
  std::stringstream sout;
  sout << "_" << g_tower.id;
  const G4String sTowerID(sout.str());

  //Processed PostionSeeds 1 from 1 1

  G4Trap* block_solid = new G4Trap(
      /*const G4String& pName*/ G4String(GetName()) + sTowerID,
      g_tower.pDz * cm,                                         // G4double pDz,
      g_tower.pTheta * rad, g_tower.pPhi * rad,                 // G4double pTheta, G4double pPhi,
      g_tower.pDy1 * cm, g_tower.pDx1 * cm, g_tower.pDx2 * cm,  // G4double pDy1, G4double pDx1, G4double pDx2,
      g_tower.pAlp1 * rad,                                      // G4double pAlp1,
      g_tower.pDy2 * cm, g_tower.pDx3 * cm, g_tower.pDx4 * cm,  // G4double pDy2, G4double pDx3, G4double pDx4,
      g_tower.pAlp2 * rad                                       // G4double pAlp2 //
  );

  G4Material* cylinder_mat = GetDetectorMaterial(get_geom_v3()->get_absorber_mat());
  assert(cylinder_mat);

  G4LogicalVolume* block_logic = new G4LogicalVolume(block_solid, cylinder_mat,
                                                     G4String(G4String(GetName()) + std::string("_Tower") + sTowerID), nullptr, nullptr,
                                                     nullptr);

  GetDisplayAction()->AddVolume(block_logic, "Block");

  // construct fibers

  if (get_geom_v3()->get_config() == SpacalGeom_t::kFullProjective_2DTaper)
  {
    int fiber_count = Construct_Fibers(g_tower, block_logic);

    if (get_geom_v3()->get_construction_verbose() >= 2)
      std::cout << "PHG4FullProjSpacalDetector::Construct_Tower::" << GetName()
                << " - constructed tower ID " << g_tower.id << " with "
                << fiber_count << " fibers using Construct_Fibers" << std::endl;
  }
  else if (get_geom_v3()->get_config() == SpacalGeom_t::kFullProjective_2DTaper_SameLengthFiberPerTower)
  {
    int fiber_count = Construct_Fibers_SameLengthFiberPerTower(g_tower,
                                                               block_logic);

    if (get_geom_v3()->get_construction_verbose() >= 2)
      std::cout << "PHG4FullProjSpacalDetector::Construct_Tower::" << GetName()
                << " - constructed tower ID " << g_tower.id << " with "
                << fiber_count
                << " fibers using Construct_Fibers_SameLengthFiberPerTower" << std::endl;
  }
  else
  {
    G4ExceptionDescription message;
    message << "can not recognize configuration type " << get_geom_v3()->get_config();

    G4Exception("PHG4FullProjSpacalDetector::Construct_Tower", "Wrong",
                FatalException, message, "");
  }

  return block_logic;
}

void PHG4FullProjSpacalDetector::Print(const std::string& /*what*/) const
{
  std::cout << "PHG4FullProjSpacalDetector::Print::" << GetName()
            << " - Print Geometry:" << std::endl;
  get_geom_v3()->Print();

  return;
}
