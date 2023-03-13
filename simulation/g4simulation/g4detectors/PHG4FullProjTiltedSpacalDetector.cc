// $$Id: PHG4FullProjTiltedSpacalDetector.cc,v 1.3 2015/02/10 15:39:26 pinkenbu Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.3 $$
 * \date $$Date: 2015/02/10 15:39:26 $$
 */
#include "PHG4FullProjTiltedSpacalDetector.h"

#include "PHG4SpacalDisplayAction.h"

#include <g4gdml/PHG4GDMLConfig.hh>

#include <phool/recoConsts.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4DisplacedSolid.hh>
#include <Geant4/G4Exception.hh>          // for G4Exception
#include <Geant4/G4ExceptionSeverity.hh>  // for FatalException
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D, G4TranslateY3D
#include <Geant4/G4Trap.hh>
#include <Geant4/G4Types.hh>  // for G4double
#include <Geant4/G4Vector3D.hh>

#include <TSystem.h>

#include <boost/foreach.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>  // for operator<<, basic_ostream
#include <limits>    // for numeric_limits
#include <map>       // for map<>::value_type, map
#include <memory>    // for allocator_traits<>::value_...
#include <numeric>   // std::accumulate
#include <sstream>
#include <string>  // std::string, std::to_string
#include <vector>  // for vector

class G4Material;
class G4VSolid;
class PHCompositeNode;

//_______________________________________________________________
//note this inactive thickness is ~1.5% of a radiation length
PHG4FullProjTiltedSpacalDetector::PHG4FullProjTiltedSpacalDetector(PHG4Subsystem* subsys, PHCompositeNode* Node,
                                                                   const std::string& dnam, PHParameters* parameters, const int lyr)
  : PHG4SpacalDetector(subsys, Node, dnam, parameters, lyr, false)
{
  assert(_geom == nullptr);

  _geom = new SpacalGeom_t();
  if (_geom == nullptr)
  {
    std::cout
        << "PHG4FullProjTiltedSpacalDetector::Constructor - Fatal Error - invalid geometry object!"
        << std::endl;
    gSystem->Exit(1);
  }
  assert(parameters);

  assert(parameters);
  get_geom_v3()->ImportParameters(*parameters);
  //  std::cout <<"PHG4FullProjTiltedSpacalDetector::Constructor -  get_geom_v3()->Print();"<<std::endl;
  //  get_geom_v3()->Print();
}

//_______________________________________________________________
void PHG4FullProjTiltedSpacalDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if (get_geom_v3()->get_construction_verbose() >= 1)
  {
    std::cout << "PHG4FullProjTiltedSpacalDetector::Construct::" << GetName()
              << " - start with PHG4SpacalDetector::Construct()." << std::endl;
  }

  PHG4SpacalDetector::ConstructMe(logicWorld);

  if (get_geom_v3()->get_construction_verbose() >= 1)
  {
    std::cout << "PHG4FullProjTiltedSpacalDetector::Construct::" << GetName()
              << " - Completed." << std::endl;
  }
}

std::pair<G4LogicalVolume*, G4Transform3D>
PHG4FullProjTiltedSpacalDetector::Construct_AzimuthalSeg()
{
  if (!(get_geom_v3()->get_azimuthal_n_sec() > 4))
  {
    std::cout << "azimuthal n sec <= 4: " << get_geom_v3()->get_azimuthal_n_sec() << std::endl;
    gSystem->Exit(1);
  }

  // basic tilt geometry
  const G4double half_chord_backend =
      get_geom_v3()->get_max_radius() * cm * tan(pi / get_geom_v3()->get_azimuthal_n_sec())  //
      + fabs(get_geom_v3()->get_thickness() * cm * 0.5 * tan(get_geom_v3()->get_azimuthal_tilt()));
  const G4double reduced_outer_radius = sqrt(pow(get_geom_v3()->get_max_radius() * cm, 2) - half_chord_backend * half_chord_backend);
  const G4double enclosure_depth = reduced_outer_radius - get_geom_v3()->get_radius() * cm;
  const G4double enclosure_center = 0.5 * (reduced_outer_radius + get_geom_v3()->get_radius() * cm);
  const G4double enclosure_half_height_half_width = enclosure_center * tan(pi / get_geom_v3()->get_azimuthal_n_sec());

  const G4double width_adj1 = tan(get_geom_v3()->get_azimuthal_tilt() - pi / get_geom_v3()->get_azimuthal_n_sec()) * enclosure_depth * 0.5;
  const G4double width_adj2 = tan(get_geom_v3()->get_azimuthal_tilt() + pi / get_geom_v3()->get_azimuthal_n_sec()) * enclosure_depth * 0.5;

  const G4double center_adj = (width_adj1 + width_adj2) * 0.5;
  const G4double center_tilt_angle = atan2(center_adj, enclosure_depth * 0.5);
  const G4double inner_half_width = enclosure_half_height_half_width + 0.5 * (width_adj1 - width_adj2);
  const G4double outter_half_width = enclosure_half_height_half_width + 0.5 * (-width_adj1 + width_adj2);

  // enclosure walls
  const G4double edge1_tilt_angle = atan2(width_adj1, enclosure_depth * 0.5);
  const G4double edge2_tilt_angle = atan2(width_adj2, enclosure_depth * 0.5);
  const G4double edge1_half_depth = sqrt(width_adj1 * width_adj1 + enclosure_depth * enclosure_depth * .25);
  const G4double edge2_half_depth = sqrt(width_adj2 * width_adj2 + enclosure_depth * enclosure_depth * .25);

  // projective center
  const G4double half_projection_ratio = 0.5 * (-width_adj1 + width_adj2) / enclosure_half_height_half_width;
  const G4double projection_center_y = enclosure_center - ((enclosure_depth * 0.5) / half_projection_ratio);
  const G4double projection_center_x = center_adj / half_projection_ratio;

  // blocks azimuthal segmentation
  const int phi_bin_in_sec = get_geom_v3()->get_max_phi_bin_in_sec();
  assert(phi_bin_in_sec >= 1);
  const G4double block_azimuth_angle = (edge2_tilt_angle - edge1_tilt_angle) / phi_bin_in_sec;
  assert(block_azimuth_angle > 0);
  if (!(fabs(block_azimuth_angle - M_PI * 2 / get_geom_v3()->get_azimuthal_n_sec() / phi_bin_in_sec) < M_PI * std::numeric_limits<G4double>::epsilon()))
  {
    std::cout << "angle/nsec out of range: " << M_PI * std::numeric_limits<G4double>::epsilon() << std::endl;
    gSystem->Exit(1);
  }
  const G4double block_edge1_half_width = enclosure_half_height_half_width - (get_geom_v3()->get_sidewall_thickness() * cm + get_geom_v3()->get_sidewall_outer_torr() * cm + 2.0 * get_geom_v3()->get_assembly_spacing() * cm) / cos(edge1_tilt_angle);
  const G4double block_edge2_half_width = enclosure_half_height_half_width - (get_geom_v3()->get_sidewall_thickness() * cm + get_geom_v3()->get_sidewall_outer_torr() * cm + 2.0 * get_geom_v3()->get_assembly_spacing() * cm) / cos(edge2_tilt_angle);
  G4double block_width_ratio = 0;
  for (int sa = 0; sa < phi_bin_in_sec; ++sa)
  {
    block_width_ratio += 1 / cos(block_azimuth_angle * (0.5 + sa) + edge1_tilt_angle);
  }
  const G4double block_half_height_width = (block_edge1_half_width + block_edge2_half_width) / block_width_ratio;
  assert(block_half_height_width > 0);

  // write out the azimuthal block geometry
  // block azimuth geometry records
  struct block_azimuth_geom
  {
    G4double angle;
    G4double projection_center_y;
    G4double projection_center_x;
    G4double projection_length;
  };
  std::vector<block_azimuth_geom> block_azimuth_geoms(phi_bin_in_sec,
                                                      block_azimuth_geom{
                                                          std::numeric_limits<double>::signaling_NaN(),
                                                          std::numeric_limits<double>::signaling_NaN(),
                                                          std::numeric_limits<double>::signaling_NaN(),
                                                          std::numeric_limits<double>::signaling_NaN()});  // [phi-bin in sector] -> azimuth geometry
  G4double block_x_edge1 = block_edge1_half_width;
  for (int sa = 0; sa < phi_bin_in_sec; ++sa)
  {
    block_azimuth_geom& geom = block_azimuth_geoms[sa];

    geom.angle = block_azimuth_angle * (0.5 + sa) + edge1_tilt_angle;
    const G4double block_x_size = block_half_height_width / cos(geom.angle);
    assert(block_x_size > 0);
    const G4double x_center = block_x_edge1 - 0.5 * block_x_size;

    // projection center per block
    geom.projection_length = block_half_height_width / 2. / tan(block_azimuth_angle / 2.);
    assert(geom.projection_length > 0);
    geom.projection_center_y = enclosure_center - geom.projection_length * cos(geom.angle);
    geom.projection_center_x = x_center + geom.projection_length * sin(geom.angle);

    // next step
    block_x_edge1 -= block_x_size;
  }

  //write out the azimuthal block divider's geometry
  struct block_divider_azimuth_geom
  {
    G4double angle;  //! rotation angle
    G4double projection_center_y;
    G4double projection_center_x;
    G4double thickness;            // thickness in the approximate azimuth direction
    G4double radial_displacement;  //! displacement along the width direction, which is the radial direction if tilt = 0
    G4double width;                //! wdith along the approximate radial direction
  };
  assert(phi_bin_in_sec >= 1);
  std::vector<block_divider_azimuth_geom> divider_azimuth_geoms(phi_bin_in_sec - 1,
                                                                block_divider_azimuth_geom{
                                                                    std::numeric_limits<double>::signaling_NaN(),
                                                                    std::numeric_limits<double>::signaling_NaN(),
                                                                    std::numeric_limits<double>::signaling_NaN(),
                                                                    std::numeric_limits<double>::signaling_NaN(),
                                                                    std::numeric_limits<double>::signaling_NaN(),
                                                                    std::numeric_limits<double>::signaling_NaN()});

  if (get_geom_v3()->get_sidewall_thickness() > 0)
  {
    for (int sa = 0; sa < phi_bin_in_sec - 1; ++sa)
    {
      block_divider_azimuth_geom& geom = divider_azimuth_geoms[sa];

      geom.angle = 0.5 * (block_azimuth_geoms[sa].angle + block_azimuth_geoms[sa + 1].angle);
      geom.projection_center_y = 0.5 * (block_azimuth_geoms[sa].projection_center_y + block_azimuth_geoms[sa + 1].projection_center_y);
      geom.projection_center_x = 0.5 * (block_azimuth_geoms[sa].projection_center_x + block_azimuth_geoms[sa + 1].projection_center_x);
      geom.radial_displacement = 0.5 * (block_azimuth_geoms[sa].projection_length + block_azimuth_geoms[sa + 1].projection_length);

      geom.thickness = 2.0 * get_geom_v3()->get_assembly_spacing() * cm * cos(block_azimuth_angle / 2.) - 2 * um;
      geom.width = get_geom_v3()->get_divider_width() * cm;
    }
  }

  if (fabs(block_x_edge1 - (-block_edge2_half_width)) > get_geom_v3()->get_assembly_spacing() * cm)
  {
    std::cout << "PHG4FullProjTiltedSpacalDetector::Construct_AzimuthalSeg - ERROR - " << std::endl
              << "\t block_x_edge1 = " << block_x_edge1 << std::endl
              << "\t block_edge2_half_width = " << block_edge2_half_width << std::endl
              << "\t fabs(block_x_edge1 - (-block_edge2_half_width)) = " << fabs(block_x_edge1 - (-block_edge2_half_width)) << std::endl
              << "\t get_geom_v3()->get_assembly_spacing() * cm = " << get_geom_v3()->get_assembly_spacing() * cm << std::endl;
  }
  if (!(fabs(block_x_edge1 - (-block_edge2_half_width)) < get_geom_v3()->get_assembly_spacing() * cm))  // closure check
  {
    std::cout << "closure check failed: " << fabs(block_x_edge1 - (-block_edge2_half_width)) << std::endl;
    gSystem->Exit(1);
  }

  if (Verbosity())
  {
    std::cout << "PHG4FullProjTiltedSpacalDetector::Construct_AzimuthalSeg - " << std::endl
              << "\t edge1_tilt_angle = " << edge1_tilt_angle << std::endl
              << "\t edge2_tilt_angle = " << edge2_tilt_angle << std::endl
              << "\t projection_center_y = " << projection_center_y << std::endl
              << "\t projection_center_x = " << projection_center_x << std::endl
              << "\t block_azimuth_angle = " << block_azimuth_angle << std::endl
              << "\t block_edge1_half_width = " << block_edge1_half_width << std::endl
              << "\t block_edge2_half_width = " << block_edge2_half_width << std::endl
              << "\t block_width_ratio = " << block_width_ratio << std::endl
              << "\t block_half_height_width = " << block_half_height_width << std::endl;

    for (int sa = 0; sa < phi_bin_in_sec; ++sa)
    {
      std::cout << "\t block[" << sa << "].angle = " << block_azimuth_geoms[sa].angle << std::endl;
      std::cout << "\t block[" << sa << "].projection_center_y = " << block_azimuth_geoms[sa].projection_center_y << std::endl;
      std::cout << "\t block[" << sa << "].projection_center_x = " << block_azimuth_geoms[sa].projection_center_x << std::endl;
    }
    for (int sa = 0; sa < phi_bin_in_sec - 1; ++sa)
    {
      std::cout << "\t divider[" << sa << "].angle = " << divider_azimuth_geoms[sa].angle << std::endl;
      std::cout << "\t divider[" << sa << "].projection_center_x = " << divider_azimuth_geoms[sa].projection_center_x << std::endl;
      std::cout << "\t divider[" << sa << "].projection_center_y = " << divider_azimuth_geoms[sa].projection_center_y << std::endl;
      std::cout << "\t divider[" << sa << "].radial_displacement = " << divider_azimuth_geoms[sa].radial_displacement << std::endl;
      std::cout << "\t divider[" << sa << "].thickness = " << divider_azimuth_geoms[sa].thickness << std::endl;
      std::cout << "\t divider[" << sa << "].width = " << divider_azimuth_geoms[sa].width << std::endl;
    }
  }

  assert(enclosure_depth > 10 * cm);

  G4VSolid* sec_solid = new G4Trap(
      /*const G4String& pName*/ G4String(GetName() + std::string("_sec_trap")),
      enclosure_depth * 0.5,                                                                              // G4double pDz,
      center_tilt_angle, halfpi,                                                                          // G4double pTheta, G4double pPhi,
      inner_half_width, get_geom_v3()->get_length() * cm / 2.0, get_geom_v3()->get_length() * cm / 2.0,   // G4double pDy1, G4double pDx1, G4double pDx2,
      0,                                                                                                  // G4double pAlp1,
      outter_half_width, get_geom_v3()->get_length() * cm / 2.0, get_geom_v3()->get_length() * cm / 2.0,  // G4double pDy2, G4double pDx3, G4double pDx4,
      0                                                                                                   // G4double pAlp2 //
  );
  G4Transform3D sec_solid_transform = G4TranslateY3D(enclosure_center) * G4RotateY3D(halfpi) * G4RotateX3D(-halfpi);
  G4VSolid* sec_solid_place = new G4DisplacedSolid(G4String(GetName() + std::string("_sec")), sec_solid, sec_solid_transform);

  recoConsts* rc = recoConsts::instance();
  G4Material* cylinder_mat = GetDetectorMaterial(rc->get_StringFlag("WorldMaterial"));
  assert(cylinder_mat);

  G4LogicalVolume* sec_logic = new G4LogicalVolume(sec_solid_place, cylinder_mat,
                                                   G4String(G4String(GetName() + std::string("_sec"))), nullptr, nullptr, nullptr);

  GetDisplayAction()->AddVolume(sec_logic, "Sector");

  // construct walls

  G4Material* wall_mat = GetDetectorMaterial(get_geom_v3()->get_sidewall_mat());
  assert(wall_mat);

  if (get_geom_v3()->get_sidewall_thickness() > 0)
  {
    // end walls
    if (get_geom_v3()->get_construction_verbose() >= 1)
    {
      std::cout << "PHG4FullProjTiltedSpacalDetector::Construct_AzimuthalSeg::" << GetName()
                << " - construct end walls." << std::endl;
    }
    //        G4Tubs* wall_solid = new G4Tubs(G4String(GetName() + std::string("_EndWall")),
    //            get_geom_v3()->get_radius() * cm + get_geom_v3()->get_sidewall_outer_torr() * cm,
    //            get_geom_v3()->get_max_radius() * cm - get_geom_v3()->get_sidewall_outer_torr() * cm,
    //            get_geom_v3()->get_sidewall_thickness() * cm / 2.0,
    //            halfpi - pi / get_geom_v3()->get_azimuthal_n_sec(),
    //            twopi / get_geom_v3()->get_azimuthal_n_sec());
    const G4double side_wall_half_thickness = get_geom_v3()->get_sidewall_thickness() * cm / 2.0;
    G4VSolid* wall_solid = new G4Trap(G4String(GetName() + std::string("_EndWall_trap")),
                                      enclosure_depth * 0.5,                                                  // G4double pDz,
                                      center_tilt_angle, halfpi,                                              // G4double pTheta, G4double pPhi,
                                      inner_half_width, side_wall_half_thickness, side_wall_half_thickness,   // G4double pDy1, G4double pDx1, G4double pDx2,
                                      0,                                                                      // G4double pAlp1,
                                      outter_half_width, side_wall_half_thickness, side_wall_half_thickness,  // G4double pDy2, G4double pDx3, G4double pDx4,
                                      0                                                                       // G4double pAlp2 //
    );
    G4VSolid* wall_solid_place = new G4DisplacedSolid(G4String(GetName() + std::string("_EndWall")), wall_solid, sec_solid_transform);

    G4LogicalVolume* wall_logic = new G4LogicalVolume(wall_solid_place, wall_mat,
                                                      G4String(G4String(GetName() + std::string("_EndWall"))), nullptr, nullptr,
                                                      nullptr);
    GetDisplayAction()->AddVolume(wall_logic, "Wall");

    using z_locations_t = std::map<int, double>;
    z_locations_t z_locations;
    z_locations[000] = get_geom_v3()->get_sidewall_thickness() * cm / 2.0 + get_geom_v3()->get_assembly_spacing() * cm;
    z_locations[001] = get_geom_v3()->get_length() * cm / 2.0 - (get_geom_v3()->get_sidewall_thickness() * cm / 2.0 + get_geom_v3()->get_assembly_spacing() * cm);
    z_locations[100] = -(get_geom_v3()->get_sidewall_thickness() * cm / 2.0 + get_geom_v3()->get_assembly_spacing() * cm);
    z_locations[101] = -(get_geom_v3()->get_length() * cm / 2.0 - (get_geom_v3()->get_sidewall_thickness() * cm / 2.0 + get_geom_v3()->get_assembly_spacing() * cm));

    BOOST_FOREACH (z_locations_t::value_type& val, z_locations)
    {
      if (get_geom_v3()->get_construction_verbose() >= 2)
      {
        std::cout << "PHG4FullProjTiltedSpacalDetector::Construct_AzimuthalSeg::"
                  << GetName() << " - constructed End Wall ID " << val.first
                  << " @ Z = " << val.second << std::endl;
      }
      G4Transform3D wall_trans = G4TranslateZ3D(val.second);

      G4PVPlacement* wall_phys = new G4PVPlacement(wall_trans, wall_logic,
                                                   G4String(GetName()) + G4String("_EndWall_") + std::to_string(val.first), sec_logic,
                                                   false, val.first, OverlapCheck());

      calo_vol[wall_phys] = val.first;
      assert(gdml_config);
      gdml_config->exclude_physical_vol(wall_phys);
    }
  }
  //
  if (get_geom_v3()->get_sidewall_thickness() > 0)
  {
    // side walls
    if (get_geom_v3()->get_construction_verbose() >= 1)
    {
      std::cout << "PHG4FullProjTiltedSpacalDetector::Construct_AzimuthalSeg::" << GetName()
                << " - construct side walls." << std::endl;
    }

    using sign_t = std::map<int, std::pair<int, int>>;
    sign_t signs;
    signs[100] = std::make_pair(+1, +1);
    signs[101] = std::make_pair(+1, -1);
    signs[200] = std::make_pair(-1, +1);
    signs[201] = std::make_pair(-1, -1);

    BOOST_FOREACH (sign_t::value_type& val, signs)
    {
      const int sign_z = val.second.first;
      const int sign_azimuth = val.second.second;

      if (get_geom_v3()->get_construction_verbose() >= 2)
      {
        std::cout << "PHG4FullProjTiltedSpacalDetector::Construct_AzimuthalSeg::"
                  << GetName() << " - constructed Side Wall ID " << val.first
                  << " with"
                  << " Shift X = "
                  << sign_azimuth * (get_geom_v3()->get_sidewall_thickness() * cm / 2.0 + get_geom_v3()->get_sidewall_outer_torr() * cm)
                  << " Rotation Z = "
                  << sign_azimuth * pi / get_geom_v3()->get_azimuthal_n_sec()
                  << " Shift Z = " << sign_z * (get_geom_v3()->get_length() * cm / 4)
                  << std::endl;
      }
      const G4double azimuth_roate = sign_azimuth > 0 ? edge1_tilt_angle : edge2_tilt_angle;
      const G4double edge_half_depth = -get_geom_v3()->get_sidewall_thickness() * cm - get_geom_v3()->get_sidewall_outer_torr() * cm + (sign_azimuth > 0 ? edge1_half_depth : edge2_half_depth);

      G4Box* wall_solid = new G4Box(G4String(GetName() + G4String("_SideWall_") + std::to_string(val.first)),
                                    get_geom_v3()->get_sidewall_thickness() * cm / 2.0,
                                    edge_half_depth,
                                    (get_geom_v3()->get_length() / 2. - 2 * (get_geom_v3()->get_sidewall_thickness() + 2. * get_geom_v3()->get_assembly_spacing())) * cm * .5);

      G4LogicalVolume* wall_logic = new G4LogicalVolume(wall_solid, wall_mat,
                                                        G4String(G4String(GetName() + G4String("_SideWall_") + std::to_string(val.first))), nullptr, nullptr,
                                                        nullptr);
      GetDisplayAction()->AddVolume(wall_logic, "Wall");

      const G4Transform3D wall_trans =
          G4TranslateZ3D(sign_z * (get_geom_v3()->get_length() * cm / 4)) *
          G4TranslateY3D(enclosure_center) *
          G4TranslateX3D(sign_azimuth * enclosure_half_height_half_width) *
          G4RotateZ3D(azimuth_roate) *
          G4TranslateX3D(-sign_azimuth * (get_geom_v3()->get_sidewall_thickness() * cm / 2.0 + get_geom_v3()->get_sidewall_outer_torr() * cm));

      G4PVPlacement* wall_phys = new G4PVPlacement(wall_trans, wall_logic,
                                                   G4String(GetName()) + G4String("_SideWall_") + std::to_string(val.first), sec_logic,
                                                   false, val.first, OverlapCheck());

      calo_vol[wall_phys] = val.first;
      assert(gdml_config);
      gdml_config->exclude_physical_vol(wall_phys);
    }
  }

  // construct dividers
  if (get_geom_v3()->get_divider_width() > 0)
  {
    if (get_geom_v3()->get_construction_verbose() >= 1)
    {
      std::cout << "PHG4FullProjTiltedSpacalDetector::Construct_AzimuthalSeg::" << GetName()
                << " - construct dividers" << std::endl;
    }

    G4Material* divider_mat = GetDetectorMaterial(get_geom_v3()->get_divider_mat());
    assert(divider_mat);

    int ID = 300;
    for (const auto& geom : divider_azimuth_geoms)
    {
      G4Box* divider_solid = new G4Box(G4String(GetName() + G4String("_Divider_") + std::to_string(ID)),
                                       geom.thickness / 2.0,
                                       geom.width / 2.,
                                       (get_geom_v3()->get_length() / 2. - 2 * (get_geom_v3()->get_sidewall_thickness() + 2. * get_geom_v3()->get_assembly_spacing())) * cm * .5);

      G4LogicalVolume* wall_logic = new G4LogicalVolume(divider_solid, divider_mat,
                                                        G4String(G4String(GetName() + G4String("_Divider_") + std::to_string(ID))), nullptr, nullptr,
                                                        nullptr);
      GetDisplayAction()->AddVolume(wall_logic, "Divider");

      for (int sign_z = -1; sign_z <= 1; sign_z += 2)
      {
        G4Transform3D wall_trans =
            G4TranslateX3D(geom.projection_center_x) *
            G4TranslateY3D(geom.projection_center_y) *
            G4RotateZ3D(geom.angle) *
            G4TranslateY3D(geom.radial_displacement) *
            G4TranslateZ3D(sign_z * (get_geom_v3()->get_length() * cm / 4));

        G4PVPlacement* wall_phys = new G4PVPlacement(wall_trans, wall_logic,
                                                     G4String(GetName()) + G4String("_Divider_") + std::to_string(ID), sec_logic,
                                                     false, ID, OverlapCheck());

        calo_vol[wall_phys] = ID;
        assert(gdml_config);
        gdml_config->exclude_physical_vol(wall_phys);

        if (Verbosity())
        {
          std::cout << "PHG4FullProjTiltedSpacalDetector::Construct_AzimuthalSeg - placing divider " << wall_phys->GetName() << " copy ID " << ID << std::endl;
        }

        ++ID;
      }
    }  //    for (const auto & geom : divider_azimuth_geoms)
  }

  //  // construct towers
  //
  BOOST_FOREACH (const SpacalGeom_t::tower_map_t::value_type& val, get_geom_v3()->get_sector_tower_map())
  {
    SpacalGeom_t::geom_tower g_tower = val.second;

    const int tower_id = g_tower.id;
    const int tower_phi_id_in_sec = tower_id % 10;
    assert(tower_phi_id_in_sec >= 0);
    assert(tower_phi_id_in_sec < phi_bin_in_sec);

    const auto& block_azimuth_geom = block_azimuth_geoms.at(tower_phi_id_in_sec);

    G4LogicalVolume* LV_tower = Construct_Tower(g_tower);

    G4Transform3D block_trans =
        G4TranslateX3D(block_azimuth_geom.projection_center_x) *
        G4TranslateY3D(block_azimuth_geom.projection_center_y) *
        G4RotateZ3D(block_azimuth_geom.angle) *
        G4TranslateX3D(g_tower.centralX * cm) *
        G4TranslateY3D(g_tower.centralY * cm) *
        G4TranslateZ3D(g_tower.centralZ * cm) *
        G4RotateX3D(g_tower.pRotationAngleX * rad);

    const bool overlapcheck_block = OverlapCheck() and (get_geom_v3()->get_construction_verbose() >= 2);

    G4PVPlacement* block_phys = new G4PVPlacement(block_trans, LV_tower,
                                                  G4String(GetName()) + G4String("_Tower_") + std::to_string(g_tower.id), sec_logic, false,
                                                  g_tower.id, overlapcheck_block);
    block_vol[block_phys] = g_tower.id;
    assert(gdml_config);
    gdml_config->exclude_physical_vol(block_phys);

    if (g_tower.LightguideHeight > 0)
    {
      // also build a light guide

      for (int ix = 0; ix < g_tower.NSubtowerX; ix++)
      //  int ix = 0;
      {
        for (int iy = 0; iy < g_tower.NSubtowerY; iy++)
        //        int iy = 0;
        {
          G4LogicalVolume* LV_lg = Construct_LightGuide(g_tower, ix,
                                                        iy);

          G4PVPlacement* lg_phys = new G4PVPlacement(block_trans, LV_lg, LV_lg->GetName(),
                                                     sec_logic, false, g_tower.id, overlapcheck_block);

          block_vol[lg_phys] = g_tower.id * 100 + ix * 10 + iy;

          assert(gdml_config);
          gdml_config->exclude_physical_vol(lg_phys);
        }
      }
    }
  }

  std::cout << "PHG4FullProjTiltedSpacalDetector::Construct_AzimuthalSeg::" << GetName()
            << " - constructed " << get_geom_v3()->get_sector_tower_map().size()
            << " unique towers" << std::endl;

  return std::make_pair(sec_logic, G4Transform3D::Identity);
}

//! Fully projective spacal with 2D tapered modules. To speed up construction, same-length fiber is used cross one tower
int PHG4FullProjTiltedSpacalDetector::Construct_Fibers_SameLengthFiberPerTower(
    const PHG4FullProjTiltedSpacalDetector::SpacalGeom_t::geom_tower& g_tower,
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
      std::cout << "PHG4FullProjTiltedSpacalDetector::Construct_Fibers_SameLengthFiberPerTower::" << GetName()
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
        << "PHG4FullProjTiltedSpacalDetector::Construct_Fibers_SameLengthFiberPerTower::"
        << GetName() << " - constructed tower ID " << g_tower.id << " with "
        << fiber_count << " fibers. Average fiber length cut = "
        << std::accumulate(fiber_cut.begin(), fiber_cut.end(), 0.0) / fiber_cut.size() << " mm" << std::endl;

  return fiber_count;
}

//! a block along z axis built with G4Trd that is slightly tapered in x dimension
int PHG4FullProjTiltedSpacalDetector::Construct_Fibers(
    const PHG4FullProjTiltedSpacalDetector::SpacalGeom_t::geom_tower& g_tower,
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
        std::cout << "PHG4FullProjTiltedSpacalDetector::Construct_Fibers::" << GetName()
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
    std::cout << "PHG4FullProjTiltedSpacalDetector::Construct_Fibers::" << GetName()
              << " - constructed tower ID " << g_tower.id << " with " << fiber_cnt
              << " fibers" << std::endl;

  return fiber_cnt;
}

//! a block along z axis built with G4Trd that is slightly tapered in x dimension
G4LogicalVolume*
PHG4FullProjTiltedSpacalDetector::Construct_Tower(
    const PHG4FullProjTiltedSpacalDetector::SpacalGeom_t::geom_tower& g_tower)
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

  if (get_geom_v3()->get_config() == SpacalGeom_t::kFullProjective_2DTaper_Tilted)
  {
    int fiber_count = Construct_Fibers(g_tower, block_logic);

    if (get_geom_v3()->get_construction_verbose() >= 2)
      std::cout << "PHG4FullProjTiltedSpacalDetector::Construct_Tower::" << GetName()
                << " - constructed tower ID " << g_tower.id << " with "
                << fiber_count << " fibers using Construct_Fibers" << std::endl;
  }
  else if (get_geom_v3()->get_config() == SpacalGeom_t::kFullProjective_2DTaper_Tilted_SameLengthFiberPerTower)
  {
    int fiber_count = Construct_Fibers_SameLengthFiberPerTower(g_tower,
                                                               block_logic);

    if (get_geom_v3()->get_construction_verbose() >= 2)
      std::cout << "PHG4FullProjTiltedSpacalDetector::Construct_Tower::" << GetName()
                << " - constructed tower ID " << g_tower.id << " with "
                << fiber_count
                << " fibers using Construct_Fibers_SameLengthFiberPerTower."
                << "V = " << block_solid->GetCubicVolume() / (cm3) << "cm3, "
                << "m = " << block_logic->GetMass() / gram << "gram, "
                << "Density = " << (block_logic->GetMass() / gram) / (block_solid->GetCubicVolume() / cm3) << "g/cm3"
                << std::endl;
  }
  else
  {
    G4ExceptionDescription message;
    message << "can not recognize configuration type " << get_geom_v3()->get_config();

    G4Exception("PHG4FullProjTiltedSpacalDetector::Construct_Tower", "Wrong",
                FatalException, message, "");
  }

  return block_logic;
}

G4LogicalVolume*
PHG4FullProjTiltedSpacalDetector::Construct_LightGuide(
    const PHG4FullProjTiltedSpacalDetector::SpacalGeom_t::geom_tower& g_tower,
    const int index_x, const int index_y)
{
  assert(_geom);
  std::stringstream sout;
  sout << "_Lightguide_" << g_tower.id << "_" << index_x << "_" << index_y;
  const G4String sTowerID(sout.str());

  assert(g_tower.LightguideHeight > 0);

  // light guide parameters in PHENIX units
  const double weight_x1 = 1 - (double) index_y / g_tower.NSubtowerY;
  const double weight_x2 = 1 - (double) (index_y + 1) / g_tower.NSubtowerY;
  const double weight_xcenter = 1 - (double) (index_y + 0.5) / g_tower.NSubtowerY;

  assert(weight_x1 >= 0 and weight_x1 <= 1);
  assert(weight_x2 >= 0 and weight_x2 <= 1);
  assert(weight_xcenter >= 0 and weight_xcenter <= 1);

  const double lg_pDx1 = (g_tower.pDx1 * weight_x1  //
                          + g_tower.pDx2 * (1 - weight_x1)) /
                         g_tower.NSubtowerX;
  const double lg_pDx2 = (g_tower.pDx1 * weight_x2  //
                          + g_tower.pDx2 * (1 - weight_x2)) /
                         g_tower.NSubtowerX;
  const double lg_pDy1 = g_tower.pDy1 / g_tower.NSubtowerY;
  const double lg_Alp1 = atan(
      (g_tower.pDx2 - g_tower.pDx1) * (-g_tower.NSubtowerX + 1. + 2 * index_x) / (double) (g_tower.NSubtowerX) / (2. * g_tower.pDy1) + tan(g_tower.pAlp1));

  const double shift_xcenter = (g_tower.pDx1 * weight_xcenter           //
                                + g_tower.pDx2 * (1 - weight_xcenter))  //
                               *                                        //
                               (-g_tower.NSubtowerX + 1. + 2 * index_x) / (double) (g_tower.NSubtowerX);
  const double shift_ycenter = g_tower.pDy1  //
                               *             //
                               (-g_tower.NSubtowerY + 1. + 2 * index_y) / (double) (g_tower.NSubtowerY);

  G4VSolid* block_solid = new G4Trap(
      /*const G4String& pName*/ G4String(GetName()) + sTowerID,
      0.5 * g_tower.LightguideHeight * cm,  // G4double pDz,
      0 * rad, 0 * rad,                     // G4double pTheta, G4double pPhi,
      g_tower.LightguideTaperRatio * lg_pDy1 * cm,
      g_tower.LightguideTaperRatio * lg_pDx1 * cm,
      g_tower.LightguideTaperRatio * lg_pDx2 * cm,  // G4double pDy1, G4double pDx1, G4double pDx2,
      lg_Alp1 * rad,                                // G4double pAlp1,
      lg_pDy1 * cm, lg_pDx1 * cm, lg_pDx2 * cm,     // G4double pDy2, G4double pDx3, G4double pDx4,
      lg_Alp1 * rad                                 // G4double pAlp2 //
  );

  block_solid = new G4DisplacedSolid(G4String(GetName() + "_displaced"),
                                     block_solid, nullptr,                                     //
                                     G4ThreeVector(                                            //
                                         tan(g_tower.pTheta * rad) * cos(g_tower.pPhi * rad),  //
                                         tan(g_tower.pTheta * rad) * sin(g_tower.pPhi * rad),  //
                                         1) *                                                  // G4ThreeVector
                                             -(g_tower.pDz) *
                                             cm                                                       //
                                         + G4ThreeVector(shift_xcenter * cm, shift_ycenter * cm, 0)   // shit in subtower direction
                                         + G4ThreeVector(0, 0, -0.5 * g_tower.LightguideHeight * cm)  //shift in the light guide height
  );

  G4Material* cylinder_mat = GetDetectorMaterial(g_tower.LightguideMaterial);
  assert(cylinder_mat);

  G4LogicalVolume* block_logic = new G4LogicalVolume(block_solid, cylinder_mat,
                                                     G4String(G4String(GetName()) + std::string("_Tower") + sTowerID), nullptr, nullptr,
                                                     nullptr);

  GetDisplayAction()->AddMaterial("LightGuide", g_tower.LightguideMaterial);
  GetDisplayAction()->AddVolume(block_logic, "LightGuide");

  return block_logic;
}

void PHG4FullProjTiltedSpacalDetector::Print(const std::string& /*what*/) const
{
  std::cout << "PHG4FullProjTiltedSpacalDetector::Print::" << GetName()
            << " - Print Geometry:" << std::endl;
  get_geom_v3()->Print();

  return;
}
