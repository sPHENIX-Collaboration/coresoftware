// $$Id: PHG4FullProjSpacalDetector.cc,v 1.3 2015/02/10 15:39:26 pinkenbu Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.3 $$
 * \date $$Date: 2015/02/10 15:39:26 $$
 */
#include "PHG4FullProjSpacalDetector.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeom_Spacalv1.h"

#include <g4main/PHG4PhenixDetector.h>
#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <fun4all/getClass.h>

#include <Geant4/G4Trd.hh>
#include <Geant4/G4Trap.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4UnionSolid.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cassert>
#include <cmath>
#include <string>     // std::string, std::to_string
#include <sstream>
#include <boost/math/special_functions/sign.hpp>
#include <boost/foreach.hpp>

using namespace std;

//_______________________________________________________________
//note this inactive thickness is ~1.5% of a radiation length
PHG4FullProjSpacalDetector::PHG4FullProjSpacalDetector(PHCompositeNode *Node,
    const std::string &dnam, SpacalGeom_t * geom, const int lyr) :
    PHG4SpacalDetector(Node, dnam,
        dynamic_cast<PHG4SpacalDetector::SpacalGeom_t *>(geom), lyr), //
    _geom(geom) //
{

  if (_geom == NULL)
    {
      cout
          << "PHG4SpacalDetector::Constructor - Fatal Error - invalid geometry object!"
          << endl;
      exit(1);
    }

  step_limits = new G4UserLimits(_geom->get_calo_step_size() * cm);

  clading_step_limits = new G4UserLimits(
      _geom->get_fiber_clading_step_size() * cm);

  fiber_core_step_limits = new G4UserLimits(
      _geom->get_fiber_core_step_size() * cm);
}

//_______________________________________________________________
void
PHG4FullProjSpacalDetector::Construct(G4LogicalVolume* logicWorld)
{
  assert(_geom);

  if (_geom->get_construction_verbose() >= 1)
    {
      cout << "PHG4FullProjSpacalDetector::Construct::" << GetName()
          << " - start with PHG4SpacalDetector::Construct()." << endl;
    }

  PHG4SpacalDetector::Construct(logicWorld);

  if (_geom->get_construction_verbose() >= 1)
    {
      cout << "PHG4FullProjSpacalDetector::Construct::" << GetName()
          << " - Completed." << endl;
    }

}

std::pair<G4LogicalVolume *, G4Transform3D>
PHG4FullProjSpacalDetector::Construct_AzimuthalSeg()
{
  assert(_geom);
  assert(_geom->get_azimuthal_n_sec()>4);

  G4Tubs* sec_solid = new G4Tubs(G4String(GetName() + string("_sec")),
      _geom->get_radius() * cm, _geom->get_max_radius() * cm,
      _geom->get_length() * cm / 2.0,
      halfpi - pi / _geom->get_azimuthal_n_sec(),
      twopi / _geom->get_azimuthal_n_sec());

  G4Material * cylinder_mat = G4Material::GetMaterial("G4_AIR");
  assert(cylinder_mat);

  G4LogicalVolume * sec_logic = new G4LogicalVolume(sec_solid, cylinder_mat,
      G4String(G4String(GetName() + string("_sec"))), 0, 0, step_limits);

  G4VisAttributes* VisAtt = new G4VisAttributes();
  VisAtt->SetColor(.1, .1, .1, .5);
  VisAtt->SetVisibility(_geom->is_virualize_fiber());
  VisAtt->SetForceSolid(false);
  sec_logic->SetVisAttributes(VisAtt);

  BOOST_FOREACH(const SpacalGeom_t::tower_map_t::value_type val, _geom->get_sector_tower_map())
    {
      const SpacalGeom_t::geom_tower & g_tower = val.second;
      G4LogicalVolume* LV_tower = Construct_Tower(g_tower);

      G4Transform3D block_trans = G4TranslateY3D(g_tower.centralY * cm)
          * G4TranslateZ3D(g_tower.centralZ * cm)
          * G4RotateX3D(g_tower.pRotationAngleX * rad);

      G4PVPlacement * block_phys = new G4PVPlacement(block_trans, LV_tower,
          G4String(GetName().c_str()) + G4String("_Tower"), sec_logic, false,
          g_tower.id, overlapcheck);
      block_vol[block_phys] = g_tower.id;

    }

  cout << "PHG4FullProjSpacalDetector::Construct_AzimuthalSeg::" << GetName()
      << " - constructed " << _geom->get_sector_tower_map().size()
      << " unique towers" << endl;
  cout << "\t" << "_geom->is_virualize_fiber() = "
      << _geom->is_virualize_fiber() << endl;

  return make_pair(sec_logic, G4Transform3D::Identity);
}

//! a block along z axis built with G4Trd that is slightly tapered in x dimension
G4LogicalVolume*
PHG4FullProjSpacalDetector::Construct_Tower(
    const PHG4FullProjSpacalDetector::SpacalGeom_t::geom_tower & g_tower)
{
  std::stringstream sout;
  sout << "_" << g_tower.id;
  const G4String sTowerID(sout.str());

  //Processed PostionSeeds 1 from 1 1

  G4Trap* block_solid = new G4Trap(
      /*const G4String& pName*/G4String(GetName().c_str()) + sTowerID,
      g_tower.pDz * cm, // G4double pDz,
      g_tower.pTheta * rad, g_tower.pPhi * rad, // G4double pTheta, G4double pPhi,
      g_tower.pDy1 * cm, g_tower.pDx1 * cm, g_tower.pDx2 * cm, // G4double pDy1, G4double pDx1, G4double pDx2,
      g_tower.pAlp1 * rad, // G4double pAlp1,
      g_tower.pDy2 * cm, g_tower.pDx3 * cm, g_tower.pDx4 * cm, // G4double pDy2, G4double pDx3, G4double pDx4,
      g_tower.pAlp2 * rad // G4double pAlp2 //
          );

  G4Material * cylinder_mat = G4Material::GetMaterial(
      _geom->get_absorber_mat());
  assert(cylinder_mat);

  G4LogicalVolume * block_logic = new G4LogicalVolume(block_solid, cylinder_mat,
      G4String(G4String(GetName()) + string("_Tower") + sTowerID), 0, 0,
      step_limits);

  G4VisAttributes* VisAtt = new G4VisAttributes();
//  PHG4Utils::SetColour(VisAtt, "W_Epoxy");
  VisAtt->SetColor(.3, .3, .3, .5);
  VisAtt->SetVisibility((!_geom->is_azimuthal_seg_visible()));
  VisAtt->SetForceWireframe((!_geom->is_azimuthal_seg_visible()));
  block_logic->SetVisAttributes(VisAtt);

  int fiber_count = 0;
//  int x_id = 0;
//  const double mid_plane_dx = _geom->get_block_width() * cm
//      * ((_geom->get_polar_taper_ratio() - 1) * 0.5 + 1) * 0.5;
//
//  const double reg_fiber_grid_distance_taper =
//      _geom->get_reg_fiber_grid_distance_taper() * cm;
//  const double reg_fiber_grid_distance_nontaper =
//      _geom->get_reg_fiber_grid_distance_nontaper() * cm;
//
//  for (double x = -mid_plane_dx + reg_fiber_grid_distance_taper / 4.;
//      x <= +mid_plane_dx; x += reg_fiber_grid_distance_taper / 2.)
//    {
//      // fiber implant
//
//      const double rotation_angle = x / mid_plane_dx
//          * atan2(
//              _geom->get_block_width() * cm
//                  * (_geom->get_polar_taper_ratio() - 1) * 0.5,
//              _geom->get_block_depth() * cm);
//
//      const double fiber_length = (_geom->get_block_depth() * cm
//          / cos(rotation_angle)) - 2 * _geom->get_fiber_outer_r() * cm;
//
//      stringstream ss;
//      ss << "_x_" << x_id;
//      G4LogicalVolume *fiber_logic = Construct_Fiber(fiber_length, ss.str());
//
//      const double dy = row_config ? reg_fiber_grid_distance_nontaper * 0.5 : 0;
//      row_config = not row_config;
//
//      for (double y = -_geom->get_block_width() * cm * 0.5
//          + reg_fiber_grid_distance_nontaper / 4.
//          + _geom->get_assembly_spacing() * cm + dy;
//          y < _geom->get_block_width() * cm * 0.5; y +=
//              reg_fiber_grid_distance_nontaper)
//        {
//          assert(
//              y -_geom->get_fiber_outer_r() * cm >-_geom->get_block_width() * cm * 0.5);
//          assert(
//              y + _geom->get_fiber_outer_r() * cm<_geom->get_block_width() * cm * 0.5);
//
//          G4Transform3D fiber_place(
//              G4TranslateY3D(y) * G4TranslateX3D(x)
//                  * G4RotateY3D(rotation_angle));
//
//          stringstream name;
//          name << GetName() << "_fiber_" << fiber_count;
//
//          G4PVPlacement * fiber_physi = new G4PVPlacement(fiber_place,
//              fiber_logic, G4String(name.str().c_str()), block_logic, false,
//              fiber_count, overlapcheck);
//          fiber_vol[fiber_physi] = fiber_count;
//
//          fiber_count++;
//        }
//
//      x_id++;
//    }

  if (_geom->get_construction_verbose())
    cout << "PHG4FullProjSpacalDetector::Construct_Block::" << GetName()
        << " - constructed tower ID " << g_tower.id << " with " << fiber_count
        << " fibers" << endl;

  return block_logic;
}

void
PHG4FullProjSpacalDetector::Print(const std::string &what) const
{
  cout << "PHG4FullProjSpacalDetector::Print::" << GetName()
      << " - Print Geometry:" << endl;
  _geom->Print();

  return;
}
