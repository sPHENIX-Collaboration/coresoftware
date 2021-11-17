#include "PHG4BbcDetector.h"

#include "PHG4BbcDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4Polyhedra.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume

#include <cmath>
#include <iostream>  // for operator<<, endl, bas...

class G4Material;
class PHCompositeNode;

PHG4BbcDetector::PHG4BbcDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *params, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4BbcDisplayAction *>(subsys->GetDisplayAction()))
  , m_Params(params)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_SupportActiveFlag(m_Params->get_int_param("supportactive"))
{
}

//_______________________________________________________________
//_______________________________________________________________
int PHG4BbcDetector::IsInBbc(G4VPhysicalVolume *volume) const
{
  G4LogicalVolume *mylogvol = volume->GetLogicalVolume();

  if (m_ActiveFlag)
  {
    if (m_PhysLogicalVolSet.find(mylogvol) != m_PhysLogicalVolSet.end())
    {
      return 1;
    }
  }
  if (m_SupportActiveFlag)
  {
    if (m_SupportLogicalVolSet.find(mylogvol) != m_SupportLogicalVolSet.end())
    {
      return -2;
    }
  }
  return 0;
}

void PHG4BbcDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  //std::cout << PHWHERE << " Constructing BBC" << std::endl;

  // Logical Volume Info for a BBC PMT
  G4Material *Quartz = GetDetectorMaterial("G4_SILICON_DIOXIDE");

  const G4double z_bbcq[] = {-1.5*cm, 1.5*cm};
  const G4double rInner_bbcq[] = {0., 0.};
  const G4double rOuter_bbcq[] = {1.27*cm, 1.27*cm};
  G4Polyhedra *bbcq = new G4Polyhedra("bbcq", 0., 2.0 * M_PI, 6, 2, z_bbcq, rInner_bbcq, rOuter_bbcq);  // bbc quartz radiator

  G4LogicalVolume *bbcq_lv = new G4LogicalVolume(bbcq, Quartz, "Bbc_quartz");

  m_PhysLogicalVolSet.insert(bbcq_lv);
  GetDisplayAction()->AddVolume(bbcq_lv, "Bbc_quartz");

  //m_bbcz = m_Params->get_double_param("z") * cm;
  m_bbcz = 250.0 * cm;  // The midpoint of the quartz is at 250 cm

  // Locations of the 64 PMT tubes in an arm.
  // Should probably move this to a geometry object...
  float xpos[] = {
      2.45951,
      -2.45951,
      4.91902,
      0,
      -4.91902,
      7.37854,
      2.45951,
      -2.45951,
      -7.37854,
      9.83805,
      4.91902,
      0,
      -4.91902,
      -9.83805,
      7.37854,
      2.45951,
      -2.45951,
      -7.37854,
      9.83805,
      4.91902,
      -4.91902,
      -9.83805,
      12.2976,
      7.37854,
      -7.37854,
      -12.2976,
      9.83805,
      -9.83805,
      12.2976,
      7.37854,
      -7.37854,
      -12.2976,
      12.2976,
      7.37854,
      -7.37854,
      -12.2976,
      9.83805,
      -9.83805,
      12.2976,
      7.37854,
      -7.37854,
      -12.2976,
      9.83805,
      4.91902,
      -4.91902,
      -9.83805,
      7.37854,
      2.45951,
      -2.45951,
      -7.37854,
      9.83805,
      4.91902,
      0,
      -4.91902,
      -9.83805,
      7.37854,
      2.45951,
      -2.45951,
      -7.37854,
      4.91902,
      0,
      -4.91902,
      2.45951,
      -2.45951,
  };

  float ypos[] = {
      12.78,
      12.78,
      11.36,
      11.36,
      11.36,
      9.94,
      9.94,
      9.94,
      9.94,
      8.52,
      8.52,
      8.52,
      8.52,
      8.52,
      7.1,
      7.1,
      7.1,
      7.1,
      5.68,
      5.68,
      5.68,
      5.68,
      4.26,
      4.26,
      4.26,
      4.26,
      2.84,
      2.84,
      1.42,
      1.42,
      1.42,
      1.42,
      -1.42,
      -1.42,
      -1.42,
      -1.42,
      -2.84,
      -2.84,
      -4.26,
      -4.26,
      -4.26,
      -4.26,
      -5.68,
      -5.68,
      -5.68,
      -5.68,
      -7.1,
      -7.1,
      -7.1,
      -7.1,
      -8.52,
      -8.52,
      -8.52,
      -8.52,
      -8.52,
      -9.94,
      -9.94,
      -9.94,
      -9.94,
      -11.36,
      -11.36,
      -11.36,
      -12.78,
      -12.78,
  };
  const int NPMT = 64;  // No. PMTs per arm
  for (int iarm = 0; iarm < 2; iarm++)
  {
    float side = 1.0;
    if (iarm == 0) side = -1.0;

    // Add BBC PMT's
    for (int itube = 0; itube < NPMT; itube++)
    {
      // Quartz Cerenkov Radiators (active)
      new G4PVPlacement(0, G4ThreeVector(side * xpos[itube] * cm, ypos[itube] * cm, side * m_bbcz),
                        bbcq_lv, "BBC", logicWorld, false, iarm * NPMT + itube, OverlapCheck());

      // PMT bodies (inactive)
    }
  }

  // Now Build the BBC Housing

  // BBC Outer and Inner Cylindrical Shells
<<<<<<< HEAD
  G4Material *Alum = G4Material::GetMaterial("G4_Al");
  G4Tubs *bbc_outer_shell = new G4Tubs("bbc_outer_shell", 14.9*cm, 15*cm, 12.2*cm, 0, 2*M_PI);
  G4LogicalVolume *bbc_outer_shell_lv = new G4LogicalVolume(bbc_outer_shell, Alum, G4String("Bbc_Outer_Shell"));
  G4VisAttributes *bbc_outer_shell_VisAtt = new G4VisAttributes();
  bbc_outer_shell_VisAtt->SetVisibility(true);
  bbc_outer_shell_VisAtt->SetForceSolid(true);
  bbc_outer_shell_VisAtt->SetColour(G4Colour::Gray());
  bbc_outer_shell_lv->SetVisAttributes(bbcqVisAtt);
  G4Tubs *bbc_inner_shell = new G4Tubs("bbc_inner_shell",5.0*cm, 5.5*cm, 12.2*cm, 0, 2*M_PI);
  G4LogicalVolume *bbc_inner_shell_lv = new G4LogicalVolume(bbc_outer_shell, Alum, G4String("Bbc_Inner_Shell"));
  G4VisAttributes *bbc_inner_shell_VisAtt = new G4VisAttributes();
  bbc_inner_shell_VisAtt->SetVisibility(true);
  bbc_inner_shell_VisAtt->SetForceSolid(true);
  bbc_inner_shell_VisAtt->SetColour(G4Colour::Gray());
  bbc_inner_shell_lv->SetVisAttributes(bbcqVisAtt);

  G4VPhysicalVolume *outer_shell_vol[2] = {0};
  G4VPhysicalVolume *inner_shell_vol[2] = {0};

  // Place South Shells
  outer_shell_vol[0] = new G4PVPlacement(0, G4ThreeVector(0, 0, (-250+2-12.2)*cm),
          bbc_outer_shell_lv, "BBC_OUTER_SHELL", logicWorld, false, 0, OverlapCheck());
  inner_shell_vol[0] = new G4PVPlacement(0, G4ThreeVector(0, 0, (-250+2-12.2)*cm),
          bbc_inner_shell_lv, "BBC_INNER_SHELL", logicWorld, false, 0, OverlapCheck());

  // Place North Shells
  outer_shell_vol[1] = new G4PVPlacement(0, G4ThreeVector(0, 0, (250-2+12.2)*cm),
          bbc_outer_shell_lv, "BBC_OUTER_SHELL", logicWorld, false, 1, OverlapCheck());
  inner_shell_vol[1] = new G4PVPlacement(0, G4ThreeVector(0, 0, (250-2+12.2)*cm),
          bbc_inner_shell_lv, "BBC_INNER_SHELL", logicWorld, false, 0, OverlapCheck());

  outer_shell_vol[0]->GetName();
  outer_shell_vol[1]->GetName();
  inner_shell_vol[0]->GetName();
  inner_shell_vol[1]->GetName();


  // BBC Front and Back Plates
  G4Tubs *bbc_plate = new G4Tubs("bbc_fplate", 5*cm, 15*cm, 0.5*cm, 0, 2*M_PI);
  G4LogicalVolume *bbc_plate_lv = new G4LogicalVolume(bbc_plate, Alum, G4String("Bbc_Cover_Plates"));
  G4VisAttributes *bbc_plate_VisAtt = new G4VisAttributes();
  bbc_plate_VisAtt->SetVisibility(true);
  bbc_plate_VisAtt->SetForceSolid(true);
  bbc_plate_VisAtt->SetColour(G4Colour::Gray());
  bbc_plate_lv->SetVisAttributes(bbcqVisAtt);

  G4VPhysicalVolume *fplate_vol[2] = {0};   // Front Plates
  G4VPhysicalVolume *bplate_vol[2] = {0};   // Back Plates

  // Place South Plates
  fplate_vol[0] = new G4PVPlacement(0, G4ThreeVector(0, 0, (-250+2.5)*cm),
          bbc_plate_lv, "BBC_FPLATE", logicWorld, false, 0, OverlapCheck());
  bplate_vol[0] = new G4PVPlacement(0, G4ThreeVector(0, 0, (-250+2.5-24.0)*cm),
          bbc_plate_lv, "BBC_BPLATE", logicWorld, false, 0, OverlapCheck());

  // Place North Plates
  fplate_vol[1] = new G4PVPlacement(0, G4ThreeVector(0, 0, (250-2.5)*cm),
          bbc_plate_lv, "BBC_FPLATE", logicWorld, false, 1, OverlapCheck());
  bplate_vol[1] = new G4PVPlacement(0, G4ThreeVector(0, 0, (250-2.5+24.0)*cm),
          bbc_plate_lv, "BBC_BPLATE", logicWorld, false, 0, OverlapCheck());

  fplate_vol[0]->GetName();
  fplate_vol[1]->GetName();
  bplate_vol[0]->GetName();
  bplate_vol[1]->GetName();

  return;
}

void PHG4BbcDetector::Print(const std::string &what) const
{
  std::cout << "Bbc Detector:" << std::endl;
  if (what == "ALL" || what == "VOLUME")
  {
    std::cout << "Version 0.1" << std::endl;
  }
  return;
}
