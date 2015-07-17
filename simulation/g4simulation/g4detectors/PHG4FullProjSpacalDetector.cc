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
#include <sstream>
#include <boost/math/special_functions/sign.hpp>

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

//
//  {
//    //Processed PostionSeeds 1 from 1 1
//    _cylinder_solid = new G4Trap( /*const G4String& pName*/G4String(GetName().c_str()),
//     6.751951 * cm,// G4double pDz,
//      0,0,// G4double pTheta, G4double pPhi,
//      2.247390*cm, 9.653108*cm , 9.642601*cm, // G4double pDy1, G4double pDx1, G4double pDx2,
//       0,// G4double pAlp1,
//      2.567902*cm, 10.983502*cm , 10.971494*cm,// G4double pDy2, G4double pDx3, G4double pDx4,
//      0// G4double pAlp2 //
//        );
//     block_trans = G4Translate3D(0.000000*cm, 105.088944*cm , 2.483168*cm) * G4RotateX3D(-1.547066*rad);
//
////  cylinder_solid = _cylinder_solid;
//
//      G4Material * cylinder_mat = G4Material::GetMaterial(_geom->get_absorber_mat());
//      assert(cylinder_mat);
//
//      cylinder_logic = new G4LogicalVolume(_cylinder_solid, cylinder_mat,
//          G4String(GetName().c_str()), 0, 0, 0);
//      G4VisAttributes* VisAtt = new G4VisAttributes();
//      PHG4Utils::SetColour(VisAtt, "W_Epoxy");
//      VisAtt->SetVisibility(true);
//      VisAtt->SetForceSolid(
//          (not _geom->is_virualize_fiber())
//              and (not _geom->is_azimuthal_seg_visible()));
//      cylinder_logic->SetVisAttributes(VisAtt);
//
//      cylinder_physi = new G4PVPlacement(block_trans, cylinder_logic,
//          G4String(GetName().c_str()), logicWorld, false, 0, overlapcheck);
//
//    }
}

std::pair<G4LogicalVolume *, G4Transform3D>
PHG4FullProjSpacalDetector::Construct_AzimuthalSeg()
{

  assert(_geom);
  assert(_geom->get_azimuthal_n_sec()>4);

  const double sec_azimuthal_width = _geom->get_sec_azimuthal_width() * cm;

  const double sec_depth = _geom->get_sec_depth() * cm;

  G4Box* sec_solid = new G4Box(G4String(GetName() + string("_sec")),
      sec_depth * .5, sec_azimuthal_width * .5, _geom->get_length() * cm * .5);

  const int sign_tilt = boost::math::sign(_geom->get_azimuthal_tilt());

  G4Transform3D sec_trans = //
      G4TranslateY3D(sign_tilt * sec_azimuthal_width * .5) //
      * G4TranslateX3D(_geom->get_radius() * cm) //
          * G4RotateZ3D(_geom->get_azimuthal_tilt()) //
          * G4TranslateY3D(-sign_tilt * sec_azimuthal_width * .5) //
          * G4TranslateX3D(sec_depth * .5);

  G4Material * cylinder_mat = G4Material::GetMaterial("G4_AIR");
  assert(cylinder_mat);

  G4LogicalVolume * sec_logic = new G4LogicalVolume(sec_solid, cylinder_mat,
      G4String(G4String(GetName() + string("_sec"))), 0, 0, step_limits);

  G4VisAttributes* VisAtt = new G4VisAttributes();
  PHG4Utils::SetColour(VisAtt, "W_Epoxy");
  VisAtt->SetVisibility(true);
  VisAtt->SetForceSolid(_geom->is_azimuthal_seg_visible());
  sec_logic->SetVisAttributes(VisAtt);

  // fill internal of sec_logic with calorimeter blocks

  G4LogicalVolume * block_logic = Construct_Block();
  assert(block_logic);

  G4Transform3D block_trans_base = G4RotateY3D(halfpi);

  double average_error = 0;
  double max_error = 0;
  double min_error = 0;
  int copy_id = 0;
  for (int sign = +1; sign >= -1; sign -= 2)
    {

      const double z_lim =
          sign > 0 ? _geom->get_zmax() * cm : _geom->get_zmin() * cm;

      double pivot_angle = sign * _geom->get_half_polar_taper_angle();
      double z_pivot_point = sign
          * (_geom->get_block_width() * cm / cos(pivot_angle)
              + _geom->get_assembly_spacing() * cm);
      while (1)
        {
          const double dz_edge = abs(
              _geom->get_block_depth() * cm * sin(pivot_angle))
              + _geom->get_block_width() * cm * 0.5
                  * (_geom->get_polar_taper_ratio() - 1);

          if (abs(z_pivot_point) + dz_edge >= abs(z_lim))
            {

              cout << "PHG4FullProjSpacalDetector::Construct_AzimuthalSeg::"
                  << GetName() << ": - "
                  << "Finished construction of blocks towards " << sign
                  << " directions. Current ID = " << copy_id << endl;

              break;
            }

          const double angle_error = halfpi
              - atan2(_geom->get_radius() * cm,
                  z_pivot_point
                      - sign * _geom->get_block_width() * cm * 0.5
                          / cos(pivot_angle)) - pivot_angle;
          average_error += sign * angle_error;
          if (sign * angle_error > max_error)
            max_error = sign * angle_error;
          if (sign * angle_error < min_error)
            min_error = sign * angle_error;

          if (_geom->get_construction_verbose() >= 2)
            {
              cout << "PHG4FullProjSpacalDetector::Construct_AzimuthalSeg::"
                  << GetName() << ": - " << "Construction of blocks " << copy_id
                  << " at " << " z_pivot_point = " << z_pivot_point
                  << " pivot_angle = " << pivot_angle << " angle_error = "
                  << angle_error << " z_lim = " << z_lim << endl;
            }

          G4Transform3D block_trans = G4TranslateZ3D(z_pivot_point) //
          * G4TranslateX3D(-_geom->get_sec_depth() * cm * 0.5) //
              * G4RotateY3D(-pivot_angle) //
              * G4TranslateZ3D(-sign * _geom->get_block_width() * cm * 0.5) //
              * G4TranslateX3D(_geom->get_block_depth() * cm * 0.5) //
              * block_trans_base;

          stringstream name;
          name << GetName() << "_block_" << copy_id;

          G4PVPlacement * block_phys = new G4PVPlacement(block_trans,
              block_logic, G4String(name.str().c_str()), sec_logic, false,
              copy_id, overlapcheck);
          block_vol[block_phys] = copy_id;

          // move to the next one
          pivot_angle += sign * _geom->get_half_polar_taper_angle() * 2;
          z_pivot_point += sign
              * (_geom->get_block_width() * cm / cos(pivot_angle)
                  + _geom->get_assembly_spacing() * cm);
          copy_id++;

          if (_geom->is_virualize_fiber())
            {
              cout << "PHG4FullProjSpacalDetector::Construct_AzimuthalSeg::"
                  << GetName() << ": - "
                  << "Only construct one block on one side for fiber illustration"
                  << endl;
              break;
            }
        }

    }

  cout << "PHG4FullProjSpacalDetector::Construct_AzimuthalSeg::" << GetName()
      << ": Averge angle error = " << average_error / copy_id << "("
      << min_error << "< Error< " << max_error << ")" << endl;
  cout << "PHG4FullProjSpacalDetector::Construct_AzimuthalSeg::" << GetName()
      << ": Done." << endl;

  return make_pair(sec_logic, sec_trans);
}

//! a block along z axis built with G4Trd that is slightly tapered in x dimension
G4LogicalVolume*
PHG4FullProjSpacalDetector::Construct_Block()
{

  G4Trd* block_solid = new G4Trd(G4String(GetName() + string("_block")),
      _geom->get_block_width() * cm * 0.5, //      G4double pdx1,
      _geom->get_block_width() * cm * _geom->get_polar_taper_ratio() * 0.5, //      G4double pdx2,
      _geom->get_block_width() * cm * 0.5, //      G4double pdy1,
      _geom->get_block_width() * cm * 0.5, //      G4double pdy2,
      _geom->get_block_depth() * cm * 0.5 //      G4double pdz
          );

  G4Material * cylinder_mat = G4Material::GetMaterial(
      _geom->get_absorber_mat());
  assert(cylinder_mat);

  G4LogicalVolume * block_logic = new G4LogicalVolume(block_solid, cylinder_mat,
      G4String(G4String(GetName() + string("_sec"))), 0, 0, step_limits);

  G4VisAttributes* VisAtt = new G4VisAttributes();
//  PHG4Utils::SetColour(VisAtt, "W_Epoxy");
  VisAtt->SetColor(.3, .3, .3, .5);

  VisAtt->SetVisibility((!_geom->is_azimuthal_seg_visible()));
//  VisAtt->SetForceSolid((!_geom->is_azimuthal_seg_visible()) and (!_geom->is_virualize_fiber()));
//  VisAtt->SetForceSolid((!_geom->is_azimuthal_seg_visible()));
  VisAtt->SetForceWireframe(true);
  block_logic->SetVisAttributes(VisAtt);

  bool row_config = true;
  int fiber_count = 0;
  int x_id = 0;
  const double mid_plane_dx = _geom->get_block_width() * cm
      * ((_geom->get_polar_taper_ratio() - 1) * 0.5 + 1) * 0.5;

  const double reg_fiber_grid_distance_taper =
      _geom->get_reg_fiber_grid_distance_taper() * cm;
  const double reg_fiber_grid_distance_nontaper =
      _geom->get_reg_fiber_grid_distance_nontaper() * cm;

  for (double x = -mid_plane_dx + reg_fiber_grid_distance_taper / 4.;
      x <= +mid_plane_dx; x += reg_fiber_grid_distance_taper / 2.)
    {
      // fiber implant

      const double rotation_angle = x / mid_plane_dx
          * atan2(
              _geom->get_block_width() * cm
                  * (_geom->get_polar_taper_ratio() - 1) * 0.5,
              _geom->get_block_depth() * cm);

      const double fiber_length = (_geom->get_block_depth() * cm
          / cos(rotation_angle)) - 2 * _geom->get_fiber_outer_r() * cm;

      stringstream ss;
      ss << "_x_" << x_id;
      G4LogicalVolume *fiber_logic = Construct_Fiber(fiber_length, ss.str());

      const double dy = row_config ? reg_fiber_grid_distance_nontaper * 0.5 : 0;
      row_config = not row_config;

      for (double y = -_geom->get_block_width() * cm * 0.5
          + reg_fiber_grid_distance_nontaper / 4.
          + _geom->get_assembly_spacing() * cm + dy;
          y < _geom->get_block_width() * cm * 0.5; y +=
              reg_fiber_grid_distance_nontaper)
        {
          assert(
              y -_geom->get_fiber_outer_r() * cm >-_geom->get_block_width() * cm * 0.5);
          assert(
              y + _geom->get_fiber_outer_r() * cm<_geom->get_block_width() * cm * 0.5);

          G4Transform3D fiber_place(
              G4TranslateY3D(y) * G4TranslateX3D(x)
                  * G4RotateY3D(rotation_angle));

          stringstream name;
          name << GetName() << "_fiber_" << fiber_count;

          G4PVPlacement * fiber_physi = new G4PVPlacement(fiber_place,
              fiber_logic, G4String(name.str().c_str()), block_logic, false,
              fiber_count, overlapcheck);
          fiber_vol[fiber_physi] = fiber_count;

          fiber_count++;
        }

      x_id++;
    }

  cout << "PHG4FullProjSpacalDetector::Construct_Block::" << GetName()
      << " - constructed " << fiber_count << " fibers in " << x_id << "rows"
      << endl;

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
