#include "PHG4CrystalCalorimeterDetector.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeomv3.h"

#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <fun4all/getClass.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4IntersectionSolid.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4Trap.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4Box.hh>

#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4Colour.hh>

#include <cmath>
#include <sstream>

using namespace std;

//static double no_overlap = 0.00015 * cm; // added safety margin against overlaps by using same boundary between volumes

//_______________________________________________________________________
PHG4CrystalCalorimeterDetector::PHG4CrystalCalorimeterDetector( PHCompositeNode *Node, const std::string &dnam ):
  PHG4Detector(Node, dnam),
  _place_in_x(0.0*mm),
  _place_in_y(0.0*mm),
  _place_in_z(-1080.0*mm),
  _rot_in_x(0.0),
  _rot_in_y(M_PI),
  _rot_in_z(0.0),
  _rMin1(22*mm),
  _rMax1(656*mm),
  _rMin2(26*mm),
  _rMax2(775*mm),
  _dZ(180*mm),
  _sPhi(0),
  _dPhi(2*M_PI),
  _crystal_front_dx(22*mm),
  _crystal_front_dy(22*mm),
  _crystal_dz(180*mm),
  _materialCrystal( "G4_PbWO4" ),
  _active(1),
  _crystallogicnameprefix("eEcalCrystal")
{

}


//_______________________________________________________________________
PHG4CrystalCalorimeterDetector::~PHG4CrystalCalorimeterDetector()
{}


//_______________________________________________________________________
int
PHG4CrystalCalorimeterDetector::IsInCrystalCalorimeter(G4VPhysicalVolume * volume) const
{
  if (volume->GetName().find(_crystallogicnameprefix) != string::npos)
    {
      return 1;
    }
  return 0;
}


//_______________________________________________________________________
void
PHG4CrystalCalorimeterDetector::Construct( G4LogicalVolume* logicWorld )
{
  /* Create the cone envelope = 'world volume' for the crystal calorimeter */
  G4Material* Air = G4Material::GetMaterial("G4_AIR");

  G4VSolid* ecal_envelope_cone = new G4Cons("eEcal_envelope_solid",
                                            _rMin1, _rMax1,
                                            _rMin2, _rMax2,
                                            _dZ/2.,
                                            _sPhi, _dPhi );

  G4LogicalVolume* ecal_envelope_log =  new G4LogicalVolume(ecal_envelope_cone, Air, G4String("eEcal_envelope"), 0, 0, 0);

  /* Define visualization attributes for envelope cone */
  G4VisAttributes* ecalVisAtt = new G4VisAttributes();
  ecalVisAtt->SetVisibility(true);
  ecalVisAtt->SetForceSolid(false);
  ecalVisAtt->SetColour(G4Colour::Magenta());
  ecal_envelope_log->SetVisAttributes(ecalVisAtt);

  /* Define rotation attributes for envelope cone */
  G4RotationMatrix ecal_rotm;
  ecal_rotm.rotateX(_rot_in_x);
  ecal_rotm.rotateY(_rot_in_y);
  ecal_rotm.rotateZ(_rot_in_z);

  /* Place envelope cone in simulation */
  new G4PVPlacement( G4Transform3D(ecal_rotm, G4ThreeVector(_place_in_x, _place_in_y, _place_in_z) ),
                     ecal_envelope_log, "CrystalCalorimeter", logicWorld, 0, false, overlapcheck);

  /* Construct crystal calorimeter within envelope */
  ConstructCrystals(ecal_envelope_log);

  return;
}


//_______________________________________________________________________
int
PHG4CrystalCalorimeterDetector::ConstructCrystals(G4LogicalVolume* ecalenvelope)
{
  /* Simple model: Fill calorimeter envelope with rectangular crystals (front face = back face)
     which are arranged in chessboard pattern in x-y and oriented parallel to z azis */
  int n_crystals_j = 58;
  int n_crystals_k = 58;

  G4double r_min = _rMin2;
  G4double r_max = _rMax1;
  G4double crystal_r = sqrt(  _crystal_front_dx * _crystal_front_dx + _crystal_front_dy * _crystal_front_dy );

  /* define center of crystal with index j=0, k=0 */
  G4double xpos_j0_k0 = 0.0*mm - ( 0.5 * n_crystals_j - 0.5 ) * _crystal_front_dx;
  G4double ypos_j0_k0 = 0.0*mm - ( 0.5 * n_crystals_k - 0.5 ) * _crystal_front_dy;

  G4Material* material_crystal = G4Material::GetMaterial(_materialCrystal.c_str());

  G4VSolid* crystal_solid = new G4Box( G4String("eEcal_crystal"),
                                       _crystal_front_dx / 2.0,
                                       _crystal_front_dy / 2.0,
                                       _crystal_dz / 2.0 );

  G4LogicalVolume *crystal_logic = new G4LogicalVolume( crystal_solid,
                                                        material_crystal,
                                                        "eEcal_crystal",
                                                        0, 0, 0);

  G4VisAttributes *visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(true);
  visattchk->SetColour(G4Colour::Cyan());
  crystal_logic->SetVisAttributes(visattchk);

  ostringstream name;

  /* Place crystal units */
  G4RotationMatrix *Rot = new G4RotationMatrix();
  Rot->rotateX(0 * rad);
  Rot->rotateY(0 * rad);
  Rot->rotateZ(0 * rad);

  for (int idx_j = 0; idx_j < n_crystals_j; idx_j++)
    {
      for (int idx_k = 0; idx_k < n_crystals_k; idx_k++)
        {
          /* Construct unique name for crystal */
          name.str("");
          name << _crystallogicnameprefix << "_j_" << idx_j << "_k_" << idx_k;

          /* Calculate center position for crystal */
          G4double xpos_crystal_jk = xpos_j0_k0 + idx_j * _crystal_front_dx;
          G4double ypos_crystal_jk = ypos_j0_k0 + idx_k * _crystal_front_dx;
          G4ThreeVector g4vec(xpos_crystal_jk, ypos_crystal_jk, 0);

          /* check if crystal extends beyond calorimeter envelope volume */
          G4double crystal_rpos = sqrt( xpos_crystal_jk * xpos_crystal_jk + ypos_crystal_jk * ypos_crystal_jk );

          G4double crystal_r_clear_max = crystal_rpos + crystal_r;
          G4double crystal_r_clear_min = crystal_rpos - crystal_r;

          if ( crystal_r_clear_min < r_min || crystal_r_clear_max > r_max )
            continue;

          /* If crystal wihtin envelope: place it */
          new G4PVPlacement( Rot, G4ThreeVector(xpos_crystal_jk , ypos_crystal_jk, 0),
                             crystal_logic,
                             name.str().c_str(),
                             ecalenvelope,
                             0, 0, overlapcheck);
        }
    }

  return 0;
}
