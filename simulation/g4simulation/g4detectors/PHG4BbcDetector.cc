#include "PHG4BbcDetector.h"

#include "PHG4BbcDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4DisplayAction.h>
#include <g4main/PHG4Subsystem.h>

#include <phool/recoConsts.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4NistManager.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4Polyhedra.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4Trap.hh>
#include <Geant4/G4Tubs.hh>

#include <Geant4/G4Types.hh>  // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>

#include <cmath>
#include <iostream>  // for operator<<, endl, bas...
#include <vector>

class PHCompositeNode;

PHG4BbcDetector::PHG4BbcDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *params, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4BbcDisplayAction *>(subsys->GetDisplayAction()))
  , m_Params(params)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_SupportActiveFlag(m_Params->get_int_param("supportactive"))
  , front_bbcz(248 * cm)
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
  // std::cout << "PHG4BbcDetector::ConstructMe()" << std::endl;

  recoConsts *rc = recoConsts::instance();

  //========================
  // Build a single BBC PMT
  //========================

  // BBC Detector Tube (BBCD) Logical Volume (contains all parts of PMT+Radiator)
  G4Material *WorldMaterial = GetDetectorMaterial(rc->get_StringFlag("WorldMaterial"));
  G4double z_bbcd[] = {-6.15 * cm, 6.15 * cm};
  G4double len_bbcd = z_bbcd[1] - z_bbcd[0];
  G4double rin_bbcd[] = {0 * cm, 0 * cm};
  G4double rout_bbcd[] = {1.4 * cm, 1.4 * cm};
  G4Polyhedra *bbcd = new G4Polyhedra("bbcd", 0., 2 * M_PI, 6, 2, z_bbcd, rin_bbcd, rout_bbcd);
  G4LogicalVolume *bbcd_lv = new G4LogicalVolume(bbcd, WorldMaterial, G4String("Bbc_tube"));

  //
  //  Place the BBCA, BBCQ, BBCP, BBCR and BBCH in to BBCD.
  //

  // BBC Attachment Plate (BBCA)
  const G4double z_bbca[] = {-0.5 * cm, -0.201 * cm, -0.2 * cm, 0.5 * cm};
  const G4double rInner_bbca[] = {0.2 * cm, 0.2 * cm, 0.2 * cm, 0.2 * cm};
  const G4double rOuter_bbca[] = {1.4 * cm, 1.4 * cm, 1.375 * cm, 1.375 * cm};
  G4Polyhedra *bbca = new G4Polyhedra("bbca", 0., 2 * M_PI, 6, 4, z_bbca, rInner_bbca, rOuter_bbca);

  G4Material *Aluminum = GetDetectorMaterial("G4_Al");

  G4LogicalVolume *bbca_lv = new G4LogicalVolume(bbca, Aluminum, G4String("Bbc_attach_plate"));
  m_SupportLogicalVolSet.insert(bbca_lv);
  GetDisplayAction()->AddVolume(bbca_lv, "Bbc_attach_plate");

  G4double xpos = 0. * cm;
  G4double ypos = 0. * cm;
  G4double len_bbca = z_bbca[3] - z_bbca[0];
  G4double zpos = z_bbcd[0] + len_bbca * 0.5;

  G4VPhysicalVolume *bbca_phys = new G4PVPlacement(nullptr, G4ThreeVector(xpos, ypos, zpos), bbca_lv, "BBCA", bbcd_lv, false, 0);
  if (!bbca_phys)  // more to prevent compiler warnings
  {
    std::cout << "placement of BBCA failed" << std::endl;
  }

  G4Material *Quartz = GetDetectorMaterial("G4_SILICON_DIOXIDE");

  // BBC Quartz Radiator
  const G4double z_bbcq[] = {-1.5 * cm, 1.5 * cm};
  const G4double rInner_bbcq[] = {0., 0.};
  const G4double rOuter_bbcq[] = {1.27 * cm, 1.27 * cm};
  G4Polyhedra *bbcq = new G4Polyhedra("bbcq", 0., 2 * M_PI, 6, 2, z_bbcq, rInner_bbcq, rOuter_bbcq);

  G4LogicalVolume *bbcq_lv = new G4LogicalVolume(bbcq, Quartz, G4String("Bbc_quartz"));
  GetDisplayAction()->AddVolume(bbcq_lv, "Bbc_quartz");

  xpos = 0. * cm;
  ypos = 0. * cm;
  G4double len_bbcq = z_bbcq[1] - z_bbcq[0];
  zpos += len_bbca * 0.5 + len_bbcq * 0.5;

  G4VPhysicalVolume *bbcq_phys = new G4PVPlacement(nullptr, G4ThreeVector(xpos, ypos, zpos), bbcq_lv, "BBCQ", bbcd_lv, false, 0);
  if (!bbcq_phys)  // more to prevent compiler warnings
  {
    std::cout << "placement of BBCQ failed" << std::endl;
  }

  // BBC PMT
  const G4double len_bbcp = 4.4 * cm;
  const G4double rInner_bbcp = 1.09 * cm;
  const G4double rOuter_bbcp = 1.29 * cm;
  G4Tubs *bbcp = new G4Tubs("bbcp", rInner_bbcp, rOuter_bbcp, len_bbcp * 0.5, 0 * deg, 360 * deg);

  G4LogicalVolume *bbcp_lv = new G4LogicalVolume(bbcp, Quartz, G4String("Bbc_PMT"));
  GetDisplayAction()->AddVolume(bbcp_lv, "Bbc_PMT");

  xpos = 0. * cm;
  ypos = 0. * cm;
  zpos += len_bbcq * 0.5 + len_bbcp * 0.5;

  G4VPhysicalVolume *bbcp_phys = new G4PVPlacement(nullptr, G4ThreeVector(xpos, ypos, zpos), bbcp_lv, "BBCP", bbcd_lv, false, 0);
  if (!bbcp_phys)  // more to prevent compiler warnings
  {
    std::cout << "placement of BBCP failed" << std::endl;
  }

  // BBC Breeder Module
  const G4double len_bbcr = 3.9 * cm;
  const G4double rInner_bbcr = 0.0 * cm;
  const G4double rOuter_bbcr = 1.2 * cm;
  G4Tubs *bbcr = new G4Tubs("bbcr", rInner_bbcr, rOuter_bbcr, len_bbcr * 0.5, 0 * deg, 360 * deg);

  G4int natoms;
  G4int ncomponents;
  G4double density;
  G4Material *G10 = new G4Material("BBC_G10", density = 1.700 * g / cm3, ncomponents = 4);
  G4NistManager *manager = G4NistManager::Instance();
  G10->AddElement(G4NistManager::Instance()->FindOrBuildElement("Si"), natoms = 1);
  G10->AddElement(G4NistManager::Instance()->FindOrBuildElement("O"), natoms = 2);
  G10->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), natoms = 3);
  G10->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), natoms = 3);

  G4LogicalVolume *bbcr_lv = new G4LogicalVolume(bbcr, G10, "Bbc_Breeder_Module");
  GetDisplayAction()->AddVolume(bbcr_lv, "Bbc_Breeder_Module");

  xpos = 0. * cm;
  ypos = 0. * cm;
  zpos += len_bbcp * 0.5 + len_bbcr * 0.5;

  G4PVPlacement *plcmt = new G4PVPlacement(nullptr, G4ThreeVector(xpos, ypos, zpos), bbcr_lv, "BBCR", bbcd_lv, false, 0);
  if (!plcmt)  // more to prevent compiler warnings
  {
    std::cout << "placement of BBCR failed" << std::endl;
  }

  // BBC Mu Metal Shield
  const G4double z_bbch[] = {-5.65 * cm, 5.65 * cm};
  const G4double rInner_bbch[] = {1.375 * cm, 1.375 * cm};
  const G4double rOuter_bbch[] = {1.4 * cm, 1.4 * cm};
  G4Polyhedra *bbch = new G4Polyhedra("bbch", 0., 2 * M_PI, 6, 2, z_bbch, rInner_bbch, rOuter_bbch);

  G4Material *MuMetal = manager->FindOrBuildMaterial("G4_STAINLESS-STEEL");

  G4LogicalVolume *bbch_lv = new G4LogicalVolume(bbch, MuMetal, G4String("Bbc_Shield"));

  xpos = 0. * cm;
  ypos = 0. * cm;
  G4double len_bbch = z_bbch[1] - z_bbch[0];
  zpos = z_bbcd[0] + 0.3 * cm + len_bbch * 0.5;

  G4VPhysicalVolume *bbch_phys = new G4PVPlacement(nullptr, G4ThreeVector(xpos, ypos, zpos), bbch_lv, "BBCH", bbcd_lv, false, 0);
  if (!bbch_phys)  // more to prevent compiler warnings
  {
    std::cout << "placement of BBCH failed" << std::endl;
  }

  // Locations of the 64 PMT tubes in an arm.
  // These are the x,y for the south BBC.
  // The north inverts the x coordinate (x -> -x)
  // (NB: Should probably move this to a geometry object...)
  float TubeLoc[64][2] = {
      {-12.2976, 4.26},
      {-12.2976, 1.42},
      {-9.83805, 8.52},
      {-9.83805, 5.68},
      {-9.83805, 2.84},
      {-7.37854, 9.94},
      {-7.37854, 7.1},
      {-7.37854, 4.26},
      {-7.37854, 1.42},
      {-4.91902, 11.36},
      {-4.91902, 8.52},
      {-4.91902, 5.68},
      {-2.45951, 12.78},
      {-2.45951, 9.94},
      {-2.45951, 7.1},
      {0, 11.36},
      {0, 8.52},
      {2.45951, 12.78},
      {2.45951, 9.94},
      {2.45951, 7.1},
      {4.91902, 11.36},
      {4.91902, 8.52},
      {4.91902, 5.68},
      {7.37854, 9.94},
      {7.37854, 7.1},
      {7.37854, 4.26},
      {7.37854, 1.42},
      {9.83805, 8.52},
      {9.83805, 5.68},
      {9.83805, 2.84},
      {12.2976, 4.26},
      {12.2976, 1.42},
      {12.2976, -4.26},
      {12.2976, -1.42},
      {9.83805, -8.52},
      {9.83805, -5.68},
      {9.83805, -2.84},
      {7.37854, -9.94},
      {7.37854, -7.1},
      {7.37854, -4.26},
      {7.37854, -1.42},
      {4.91902, -11.36},
      {4.91902, -8.52},
      {4.91902, -5.68},
      {2.45951, -12.78},
      {2.45951, -9.94},
      {2.45951, -7.1},
      {0, -11.36},
      {0, -8.52},
      {-2.45951, -12.78},
      {-2.45951, -9.94},
      {-2.45951, -7.1},
      {-4.91902, -11.36},
      {-4.91902, -8.52},
      {-4.91902, -5.68},
      {-7.37854, -9.94},
      {-7.37854, -7.1},
      {-7.37854, -4.26},
      {-7.37854, -1.42},
      {-9.83805, -8.52},
      {-9.83805, -5.68},
      {-9.83805, -2.84},
      {-12.2976, -4.26},
      {-12.2976, -1.42}};

  // m_bbcz = m_Params->get_double_param("z") * cm;
  m_bbcz = 250.0 * cm;  // The front face of the quartz is at 250 cm

  const float tube_zpos = m_bbcz + len_bbcd / 2.0 - len_bbca;

  const int NPMT = 64;  // No. PMTs per arm
  G4RotationMatrix *arm_rot[2];
  for (int iarm = 0; iarm < 2; iarm++)
  {
    arm_rot[iarm] = new G4RotationMatrix;
    float xside = 1.0;
    float zside = 1.0;
    if (iarm == 0)  // South Arm = 0
    {
      xside = 1.0;
      zside = -1.0;
      arm_rot[iarm]->rotateY(180 * deg);
    }
    else  // North Arm = 1
    {
      xside = -1.0;
      zside = 1.0;
    }

    // Add BBC PMT's
    for (int itube = 0; itube < NPMT; itube++)
    {
      // Full PMT Housing with Active Quartz Cerenkov Radiators
      float tube_xpos = xside * TubeLoc[itube][0] * cm;
      float tube_ypos = TubeLoc[itube][1] * cm;
      new G4PVPlacement(arm_rot[iarm], G4ThreeVector(tube_xpos, tube_ypos, zside * tube_zpos),
                        bbcd_lv, "BBCD", logicWorld, false, iarm * NPMT + itube, OverlapCheck());
    }
  }

  // Now Build the BBC Housing

  // BBC Outer and Inner Cylindrical Shells
  G4Tubs *bbc_outer_shell = new G4Tubs("bbc_outer_shell", 14.9 * cm, 15 * cm, 11.5 * cm, 0, 2 * M_PI);
  G4LogicalVolume *bbc_outer_shell_lv = new G4LogicalVolume(bbc_outer_shell, Aluminum, G4String("Bbc_Outer_Shell"));
  GetDisplayAction()->AddVolume(bbc_outer_shell_lv, "Bbc_Outer_Shell");

  G4Tubs *bbc_inner_shell = new G4Tubs("bbc_inner_shell", 5.0 * cm, 5.5 * cm, 11.5 * cm, 0, 2 * M_PI);
  G4LogicalVolume *bbc_inner_shell_lv = new G4LogicalVolume(bbc_inner_shell, Aluminum, G4String("Bbc_Inner_Shell"));
  GetDisplayAction()->AddVolume(bbc_inner_shell_lv, "Bbc_Inner_Shell");

  G4VPhysicalVolume *outer_shell_vol[2] = {nullptr};
  G4VPhysicalVolume *inner_shell_vol[2] = {nullptr};

  // Place South Shells
  outer_shell_vol[0] = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, (-250 + 1.0 - 11.5) * cm),
                                         bbc_outer_shell_lv, "BBC_OUTER_SHELL", logicWorld, false, 0, OverlapCheck());
  inner_shell_vol[0] = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, (-250 + 1.0 - 11.5) * cm),
                                         bbc_inner_shell_lv, "BBC_INNER_SHELL", logicWorld, false, 0, OverlapCheck());

  // Place North Shells
  outer_shell_vol[1] = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, (250 - 1.0 + 11.5) * cm),
                                         bbc_outer_shell_lv, "BBC_OUTER_SHELL", logicWorld, false, 1, OverlapCheck());
  inner_shell_vol[1] = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, (250 - 1.0 + 11.5) * cm),
                                         bbc_inner_shell_lv, "BBC_INNER_SHELL", logicWorld, false, 0, OverlapCheck());
  // this is more to prevent compiler warnings about unused variables
  if (!outer_shell_vol[0] || !outer_shell_vol[1] || !inner_shell_vol[0] || !inner_shell_vol[1])
  {
    std::cout << "problem placing BBC Sheels" << std::endl;
  }

  // BBC Front and Back Plates
  G4Tubs *bbc_plate = new G4Tubs("bbc_fplate", 5 * cm, 15 * cm, 0.5 * cm, 0, 2 * M_PI);
  G4LogicalVolume *bbc_plate_lv = new G4LogicalVolume(bbc_plate, Aluminum, G4String("Bbc_Cover_Plates"));
  GetDisplayAction()->AddVolume(bbc_plate_lv, "Bbc_Cover_Plates");

  G4VPhysicalVolume *fplate_vol[2] = {nullptr};  // Front Plates
  G4VPhysicalVolume *bplate_vol[2] = {nullptr};  // Back Plates

  // Place South Plates
  fplate_vol[0] = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, (-250 + 1.0 + 0.5) * cm),
                                    bbc_plate_lv, "BBC_FPLATE", logicWorld, false, 0, OverlapCheck());
  bplate_vol[0] = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, (-250 + 1.0 + 0.5 - 24.0) * cm),
                                    bbc_plate_lv, "BBC_BPLATE", logicWorld, false, 0, OverlapCheck());

  // Place North Plates
  fplate_vol[1] = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, (250 - 1.0 - 0.5) * cm),
                                    bbc_plate_lv, "BBC_FPLATE", logicWorld, false, 1, OverlapCheck());
  bplate_vol[1] = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, (250 - 1.0 - 0.5 + 24.0) * cm),
                                    bbc_plate_lv, "BBC_BPLATE", logicWorld, false, 0, OverlapCheck());

  // Place BBC Cables
  G4Material *Cu = manager->FindOrBuildMaterial("G4_Cu");
  const G4double len_cable = 120 * cm;
  const G4double r_CableConductor = 0.09525 * cm;
  G4Tubs *bbc_cablecond = new G4Tubs("bbc_cablecond", 0., r_CableConductor, len_cable * 0.5, 0 * deg, 360 * deg);

  G4LogicalVolume *bbc_cablecond_lv = new G4LogicalVolume(bbc_cablecond, Cu, G4String("Bbc_CableCond"));
  GetDisplayAction()->AddVolume(bbc_cablecond_lv, "Bbc_CableCond");

  const G4double rIn_CableShield = 0.302876 * cm;
  const G4double rOut_CableShield = 0.3175 * cm;
  G4Tubs *bbc_cableshield = new G4Tubs("bbc_cableshield", rIn_CableShield, rOut_CableShield, len_cable * 0.5, 0 * deg, 360 * deg);

  G4LogicalVolume *bbc_cableshield_lv = new G4LogicalVolume(bbc_cableshield, Cu, G4String("Bbc_CableShield"));
  GetDisplayAction()->AddVolume(bbc_cableshield_lv, "Bbc_CableShield");

  ypos = len_cable / 2 + 5 * cm;

  // For now we make this vertical, but they should be slanted toward the endcap
  G4RotationMatrix *rot_cable = new G4RotationMatrix();
  rot_cable->rotateX(90 * deg);

  int icable = 0;

  for (int iarm = 0; iarm < 2; iarm++)
  {
    float zsign = -1.0;
    if (iarm == 1)
    {
      zsign = 1.0;
    }

    for (int iring = 1; iring < 5; iring++)
    {
      float ring_radius = iring * 0.67 * cm;
      int ncables = 2 * M_PI * ring_radius / (0.635 * cm);
      double dphi = 2 * M_PI / ncables;

      // G4cout << "BBC_CABLE " << iring << "\t" << ring_radius << "\t" << ncables << "\t" << dphi*180/3.14 << G4endl;

      // place cables in ring
      for (int ic = 0; ic < ncables; ic++)
      {
        xpos = cos(dphi * ic) * ring_radius;
        zpos = sin(dphi * ic) * ring_radius + zsign * (m_bbcz + 30 * cm);
        // Place Inner Conductor
        new G4PVPlacement(rot_cable, G4ThreeVector(xpos, ypos, zpos), bbc_cablecond_lv, "BBC_Cable_Cond", logicWorld, false, icable, OverlapCheck());
        // Place Shield
        new G4PVPlacement(rot_cable, G4ThreeVector(xpos, ypos, zpos), bbc_cableshield_lv, "BBC_Cable_Shield", logicWorld, false, icable++, OverlapCheck());
      }
    }
  }

  // Now make the Supports
  ConstructSupport(logicWorld);

  // this is more to prevent compiler warnings about unused variables
  if (!fplate_vol[0] || !fplate_vol[1] || !bplate_vol[0] || !bplate_vol[1])
  {
    std::cout << "problem placing BBC Sheets" << std::endl;
  }

  // bbcq is the active detector element
  m_PhysLogicalVolSet.insert(bbcq_lv);

  //  GetDisplayAction()->AddVolume(bbc_plate_lv, "Bbc_plate");
  //  GetDisplayAction()->AddVolume(bbc_outer_shell_lv, "Bbc_oshell");
  //  GetDisplayAction()->AddVolume(bbc_inner_shell_lv, "Bbc_ishell");

  std::cout << "INSIDE BBC" << std::endl;

  return;
}

void PHG4BbcDetector::ConstructSupport(G4LogicalVolume *logicWorld)
{
  // std::cout << "PHG4BbcDetector::ConstructSupport()" << std::endl;

  G4double fractionmass;
  G4int ncomponents;
  G4double density;

  recoConsts *rc = recoConsts::instance();
  G4Material *WorldMaterial = GetDetectorMaterial(rc->get_StringFlag("WorldMaterial"));
  G4Material *Delrin = GetDetectorMaterial("G4_POLYOXYMETHYLENE");
  // G4Material *Aluminum = GetDetectorMaterial("G4_Al");

  // McMaster-Carr #8548K36, implemented as PDG G10, https://pdg.lbl.gov/2022/AtomicNuclearProperties/HTML/G10.html
  G4Material *Fiberglass = new G4Material("BBC_Fiberglass", density = 1.800 * g / cm3, ncomponents = 9);
  Fiberglass->AddElement(G4NistManager::Instance()->FindOrBuildElement("B"), fractionmass = 0.018640);
  Fiberglass->AddElement(G4NistManager::Instance()->FindOrBuildElement("Mg"), fractionmass = 0.010842);
  Fiberglass->AddElement(G4NistManager::Instance()->FindOrBuildElement("O"), fractionmass = 0.044453);
  Fiberglass->AddElement(G4NistManager::Instance()->FindOrBuildElement("Al"), fractionmass = 0.151423);
  Fiberglass->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), fractionmass = 0.081496);
  Fiberglass->AddElement(G4NistManager::Instance()->FindOrBuildElement("Ca"), fractionmass = 0.385764);
  Fiberglass->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), fractionmass = 0.275853);
  Fiberglass->AddElement(G4NistManager::Instance()->FindOrBuildElement("N"), fractionmass = 0.003583);
  Fiberglass->AddElement(G4NistManager::Instance()->FindOrBuildElement("Si"), fractionmass = 0.027945);

  // BBC Base Plates
  G4double basep_width = 35.56 * cm;
  G4double basep_height = 1.91 * cm;
  G4double basep_len = 46.99 * cm;
  G4double basep_zpos = 228.6 * cm;  // z-pos of front edge of base plate
  G4Box *bbc_base_plate = new G4Box("bbc_base_plate", basep_width / 2, basep_height / 2, basep_len / 2);
  G4LogicalVolume *bbc_base_plate_lv = new G4LogicalVolume(bbc_base_plate, Delrin, G4String("Bbc_Base_Plates"));
  GetDisplayAction()->AddVolume(bbc_base_plate_lv, "Bbc_Base_Plates");

  // Place South Base Plates
  new G4PVPlacement(nullptr, G4ThreeVector(0 * cm, -15. * cm - basep_height / 2, (-basep_zpos - basep_len / 2)),
                    bbc_base_plate_lv, "BBC_BASE_PLATE", logicWorld, false, 0, OverlapCheck());

  // Place North Base Plates
  new G4PVPlacement(nullptr, G4ThreeVector(0 * cm, -15. * cm - basep_height / 2, (basep_zpos + basep_len / 2)),
                    bbc_base_plate_lv, "BBC_BASE_PLATE", logicWorld, false, 1, OverlapCheck());

  // BBC Side Support Plates
  G4double sidesupportp_width = 1.27 * cm;
  G4double sidesupportp_height = 14.57 * cm;
  G4double sidesupportp_len = 25.00 * cm;

  G4Box *bbc_sidesupport_plate = new G4Box("bbc_sidesupport_plate", sidesupportp_width / 2, sidesupportp_height / 2, sidesupportp_len / 2);
  G4LogicalVolume *bbc_sidesupport_plate_lv = new G4LogicalVolume(bbc_sidesupport_plate, Delrin, G4String("Bbc_Sidesupport_Plates"));
  GetDisplayAction()->AddVolume(bbc_sidesupport_plate_lv, "Bbc_Sidesupport_Plates");

  // Make and Place Holes in Side Support Plates
  G4Tubs *bbc_sidesupport_hole = new G4Tubs("bbc_sidesupport_hole", 0., (7.62 / 2) * cm, sidesupportp_width / 2, 0, 2.0 * M_PI);
  G4LogicalVolume *bbc_sidesupport_hole_lv = new G4LogicalVolume(bbc_sidesupport_hole, WorldMaterial, G4String("Bbc_Sidesupport_Holes"));

  G4RotationMatrix *rot_sideholes = new G4RotationMatrix;
  rot_sideholes->rotateY(90. * deg);
  new G4PVPlacement(rot_sideholes, G4ThreeVector(0, -0.935 * cm, 11.43 * cm / 2), bbc_sidesupport_hole_lv, "BBC_SIDESUPPORT_HOLE", bbc_sidesupport_plate_lv, false, 0, OverlapCheck());
  new G4PVPlacement(rot_sideholes, G4ThreeVector(0, -0.935 * cm, -11.43 * cm / 2), bbc_sidesupport_hole_lv, "BBC_SIDESUPPORT_HOLE", bbc_sidesupport_plate_lv, false, 1, OverlapCheck());

  // Place South Side Support Plates
  new G4PVPlacement(nullptr, G4ThreeVector(-15 * cm - sidesupportp_width / 2, -15 * cm + sidesupportp_height / 2, -front_bbcz - sidesupportp_len / 2),
                    bbc_sidesupport_plate_lv, "BBC_SIDESUPPORT_PLATE", logicWorld, false, 0, OverlapCheck());
  new G4PVPlacement(nullptr, G4ThreeVector(15 * cm + sidesupportp_width / 2, -15 * cm + sidesupportp_height / 2, -front_bbcz - sidesupportp_len / 2),
                    bbc_sidesupport_plate_lv, "BBC_SIDESUPPORT_PLATE", logicWorld, false, 1, OverlapCheck());

  // Place North Side Support Plate
  new G4PVPlacement(nullptr, G4ThreeVector(-15 * cm - sidesupportp_width / 2, -15 * cm + sidesupportp_height / 2, front_bbcz + sidesupportp_len / 2),
                    bbc_sidesupport_plate_lv, "BBC_SIDESUPPORT_PLATE", logicWorld, false, 2, OverlapCheck());
  new G4PVPlacement(nullptr, G4ThreeVector(15 * cm + sidesupportp_width / 2, -15 * cm + sidesupportp_height / 2, front_bbcz + sidesupportp_len / 2),
                    bbc_sidesupport_plate_lv, "BBC_SIDESUPPORT_PLATE", logicWorld, false, 3, OverlapCheck());

  // BBC Support Post
  G4double supportp_width = 10.16 * cm;
  G4double supportp_height = 97.16 * cm;
  G4double supportp_len = 10.16 * cm;
  G4double supportp_thick = 0.64 * cm;

  G4Box *bbc_support_post_outside = new G4Box("bbc_support_post_outside", supportp_width / 2, supportp_height / 2, supportp_len / 2);
  G4Box *bbc_support_post_inside = new G4Box("bbc_support_post_inside", supportp_width / 2 - supportp_thick, supportp_height / 2 - supportp_thick, supportp_len / 2 - supportp_thick);
  G4SubtractionSolid *bbc_support_post = new G4SubtractionSolid("support_post", bbc_support_post_outside, bbc_support_post_inside);

  G4LogicalVolume *bbc_support_post_lv = new G4LogicalVolume(bbc_support_post, Fiberglass, G4String("Bbc_Support_Post"));

  GetDisplayAction()->AddVolume(bbc_support_post_lv, "Bbc_Support_Post");

  // Place South Support Post
  new G4PVPlacement(nullptr, G4ThreeVector(0, -15. * cm - basep_height - supportp_height / 2, -228. * cm - supportp_len / 2),
                    bbc_support_post_lv, "BBC_SUPPORT_POST", logicWorld, false, 0, OverlapCheck());

  // Place North Support Post
  new G4PVPlacement(nullptr, G4ThreeVector(0, -15. * cm - basep_height - supportp_height / 2, 228. * cm + supportp_len / 2),
                    bbc_support_post_lv, "BBC_SUPPORT_POST", logicWorld, false, 1, OverlapCheck());

  // BBC Support Arm
  G4double supporta_width = 148.59 * cm;
  G4double supporta_height = 10.16 * cm;
  G4double supporta_len = 10.16 * cm;
  G4double supporta_thick = 0.64 * cm;

  G4Box *bbc_support_arm_outside = new G4Box("bbc_support_arm_outside", supporta_width / 2, supporta_height / 2, supporta_len / 2);
  G4Box *bbc_support_arm_inside = new G4Box("bbc_support_arm_inside", supporta_width / 2 - supporta_thick, supporta_height / 2 - supporta_thick, supporta_len / 2 - supporta_thick);
  G4SubtractionSolid *bbc_support_arm = new G4SubtractionSolid("support_arm", bbc_support_arm_outside, bbc_support_arm_inside);

  G4LogicalVolume *bbc_support_arm_lv = new G4LogicalVolume(bbc_support_arm, Fiberglass, G4String("Bbc_Support_Arm"));

  GetDisplayAction()->AddVolume(bbc_support_arm_lv, "Bbc_Support_Arm");

  // Place South Support Arms
  new G4PVPlacement(nullptr, G4ThreeVector(-supporta_width / 2 - supportp_width / 2, -15. * cm - basep_height - 20.48 * cm - supporta_height / 2, -228.6 * cm - supporta_len / 2),
                    bbc_support_arm_lv, "BBC_SUPPORT_ARM", logicWorld, false, 0, OverlapCheck());
  new G4PVPlacement(nullptr, G4ThreeVector(supporta_width / 2 + supportp_width / 2, -15. * cm - basep_height - 20.48 * cm - supporta_height / 2, -228.6 * cm - supporta_len / 2),
                    bbc_support_arm_lv, "BBC_SUPPORT_ARM", logicWorld, false, 1, OverlapCheck());

  // Place North Support Arms
  new G4PVPlacement(nullptr, G4ThreeVector(-supporta_width / 2 - supportp_width / 2, -15. * cm - basep_height - 20.48 * cm - supporta_height / 2, 228.6 * cm + supporta_len / 2),
                    bbc_support_arm_lv, "BBC_SUPPORT_ARM", logicWorld, false, 2, OverlapCheck());
  new G4PVPlacement(nullptr, G4ThreeVector(supporta_width / 2 + supportp_width / 2, -15. * cm - basep_height - 20.48 * cm - supporta_height / 2, 228.6 * cm + supporta_len / 2),
                    bbc_support_arm_lv, "BBC_SUPPORT_ARM", logicWorld, false, 3, OverlapCheck());

  // BBC Gusset Plates (implementation is broken up into 3 parts, in Trap and 2 Boxes)
  G4double gussetp0_pz = 1.27 * cm;
  G4double gussetp0_py = (45.72 - 11.11) * cm;
  G4double gussetp0_px = 44.62 * cm;  // measured off drawing
  G4double gussetp0_pltx = 5.08 * cm;
  G4Trap *bbc_gusset_plate0 = new G4Trap("bbc_gusset_plate0", gussetp0_pz, gussetp0_py, gussetp0_px, gussetp0_pltx);
  G4LogicalVolume *bbc_gusset0_plate_lv = new G4LogicalVolume(bbc_gusset_plate0, Delrin, G4String("Bbc_Gusset0_Plates"));
  GetDisplayAction()->AddVolume(bbc_gusset0_plate_lv, "Bbc_Gusset0_Plates");

  // Place South Gusset Plates (Trapezoid part)
  G4RotationMatrix *rot_sgusset = new G4RotationMatrix;
  rot_sgusset->rotateY(90. * deg);
  rot_sgusset->rotateZ(180. * deg);

  G4double xpos = supportp_width / 2 + gussetp0_pz / 2;
  G4double ypos = -gussetp0_py / 2 - 15 * cm - basep_height;
  G4double zpos = 0.25 * (gussetp0_px + gussetp0_pltx) + 228.6 * cm + 11.11 * cm;
  new G4PVPlacement(rot_sgusset, G4ThreeVector(-xpos, ypos, -zpos), bbc_gusset0_plate_lv, "BBC_GUSSET_PLATE0", logicWorld, false, 0, OverlapCheck());
  new G4PVPlacement(rot_sgusset, G4ThreeVector(xpos, ypos, -zpos), bbc_gusset0_plate_lv, "BBC_GUSSET_PLATE0", logicWorld, false, 1, OverlapCheck());

  // Place North Gusset Plates (Trapezoid part)
  G4RotationMatrix *rot_ngusset = new G4RotationMatrix;
  rot_ngusset->rotateY(-90. * deg);
  rot_ngusset->rotateZ(180. * deg);

  new G4PVPlacement(rot_ngusset, G4ThreeVector(-xpos, ypos, zpos), bbc_gusset0_plate_lv, "BBC_GUSSET_PLATE0", logicWorld, false, 2, OverlapCheck());
  new G4PVPlacement(rot_ngusset, G4ThreeVector(xpos, ypos, zpos), bbc_gusset0_plate_lv, "BBC_GUSSET_PLATE0", logicWorld, false, 3, OverlapCheck());

  // Gusset Plate Box 1
  G4double gussetp1_x = 1.27 * cm;  // measured off drawing
  G4double gussetp1_y = 20.48 * cm;
  G4double gussetp1_z = 11.11 * cm;

  G4Box *bbc_gusset_plate1 = new G4Box("bbc_gusset_plate1", gussetp1_x / 2, gussetp1_y / 2, gussetp1_z / 2);
  G4LogicalVolume *bbc_gusset1_plate_lv = new G4LogicalVolume(bbc_gusset_plate1, Delrin, G4String("Bbc_Gusset1_Plates"));
  GetDisplayAction()->AddVolume(bbc_gusset1_plate_lv, "Bbc_Gusset1_Plates");

  // Place South Gusset Plates (Box 1)
  ypos = -15 * cm - basep_height - gussetp1_y / 2;
  zpos = 228.6 * cm + gussetp1_z / 2;
  new G4PVPlacement(nullptr, G4ThreeVector(-xpos, ypos, -zpos), bbc_gusset1_plate_lv, "BBC_GUSSET_PLATE1", logicWorld, false, 0, OverlapCheck());
  new G4PVPlacement(nullptr, G4ThreeVector(xpos, ypos, -zpos), bbc_gusset1_plate_lv, "BBC_GUSSET_PLATE1", logicWorld, false, 1, OverlapCheck());

  // Place North Gusset Plates (Box 1)
  new G4PVPlacement(nullptr, G4ThreeVector(-xpos, ypos, zpos), bbc_gusset1_plate_lv, "BBC_GUSSET_PLATE1", logicWorld, false, 2, OverlapCheck());
  new G4PVPlacement(nullptr, G4ThreeVector(xpos, ypos, zpos), bbc_gusset1_plate_lv, "BBC_GUSSET_PLATE1", logicWorld, false, 3, OverlapCheck());

  // Gusset Plate Box 2
  G4double gussetp2_x = 1.27 * cm;  // measured off drawing
  G4double gussetp2_y = (45.72 - 10.80 - 20.48) * cm;
  G4double gussetp2_z = 11.11 * cm;

  G4Box *bbc_gusset_plate2 = new G4Box("bbc_gusset_plate2", gussetp2_x / 2, gussetp2_y / 2, gussetp2_z / 2);
  G4LogicalVolume *bbc_gusset2_plate_lv = new G4LogicalVolume(bbc_gusset_plate2, Delrin, G4String("Bbc_Gusset2_Plates"));
  GetDisplayAction()->AddVolume(bbc_gusset2_plate_lv, "Bbc_Gusset2_Plates");

  // Place South Gusset Plates (Box 2)
  ypos = -15 * cm - basep_height - 20.48 * cm - 10.80 * cm - gussetp2_y / 2;
  zpos = 228.6 * cm + gussetp2_z / 2;
  new G4PVPlacement(nullptr, G4ThreeVector(-xpos, ypos, -zpos), bbc_gusset2_plate_lv, "BBC_GUSSET_PLATE2", logicWorld, false, 0, OverlapCheck());
  new G4PVPlacement(nullptr, G4ThreeVector(xpos, ypos, -zpos), bbc_gusset2_plate_lv, "BBC_GUSSET_PLATE2", logicWorld, false, 1, OverlapCheck());

  // Place North Gusset Plates (Box 2)
  new G4PVPlacement(nullptr, G4ThreeVector(-xpos, ypos, zpos), bbc_gusset2_plate_lv, "BBC_GUSSET_PLATE2", logicWorld, false, 2, OverlapCheck());
  new G4PVPlacement(nullptr, G4ThreeVector(xpos, ypos, zpos), bbc_gusset2_plate_lv, "BBC_GUSSET_PLATE2", logicWorld, false, 3, OverlapCheck());

  // BBC Splice Plates
  G4double splicep_x = 35.56 * cm;
  G4double splicep_y = 10.16 * cm;
  G4double splicep_z = 0.64 * cm;

  G4Box *bbc_splice_plate = new G4Box("bbc_splice_plate", splicep_x / 2, splicep_y / 2, splicep_z / 2);
  G4LogicalVolume *bbc_splice_plate_lv = new G4LogicalVolume(bbc_splice_plate, Delrin, G4String("Bbc_Splice_Plates"));
  GetDisplayAction()->AddVolume(bbc_splice_plate_lv, "Bbc_Splice_Plates");

  // Make and Place Holes in Splice Plates
  G4Tubs *bbc_splice_hole = new G4Tubs("bbc_splice_hole", 0., (6.35 / 2) * cm, splicep_z / 2, 0, 2.0 * M_PI);
  G4LogicalVolume *bbc_splice_hole_lv = new G4LogicalVolume(bbc_splice_hole, WorldMaterial, G4String("Bbc_Splice_Holes"));

  new G4PVPlacement(nullptr, G4ThreeVector(-12.7 * cm, 0, 0), bbc_splice_hole_lv, "BBC_SPLICE_HOLE", bbc_splice_plate_lv, false, 0, OverlapCheck());
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), bbc_splice_hole_lv, "BBC_SPLICE_HOLE", bbc_splice_plate_lv, false, 1, OverlapCheck());
  new G4PVPlacement(nullptr, G4ThreeVector(12.7 * cm, 0, 0), bbc_splice_hole_lv, "BBC_SPLICE_HOLE", bbc_splice_plate_lv, false, 2, OverlapCheck());

  // Place South Splice Plate
  ypos = -15 * cm - basep_height - 20.48 * cm - splicep_y / 2;
  zpos = 228.6 * cm + supportp_len + splicep_z / 2;
  new G4PVPlacement(nullptr, G4ThreeVector(0, ypos, -zpos), bbc_splice_plate_lv, "BBC_SPLICE_PLATE", logicWorld, false, 0, OverlapCheck());

  // Place North Splice Plate
  new G4PVPlacement(nullptr, G4ThreeVector(0, ypos, zpos), bbc_splice_plate_lv, "BBC_SPLICE_PLATE", logicWorld, false, 1, OverlapCheck());
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
