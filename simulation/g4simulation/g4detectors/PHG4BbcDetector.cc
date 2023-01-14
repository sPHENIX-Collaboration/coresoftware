#include "PHG4BbcDetector.h"

#include "PHG4BbcDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4DisplayAction.h>
#include <g4main/PHG4Subsystem.h>

#include <phool/recoConsts.h>

#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4NistManager.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4Polyhedra.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Types.hh>  // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>

#include <cmath>
#include <iostream>  // for operator<<, endl, bas...

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
  //std::cout << "PHG4BbcDetector::ConstructMe()" << std::endl;

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

  //m_bbcz = m_Params->get_double_param("z") * cm;
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

      //G4cout << "BBC_CABLE " << iring << "\t" << ring_radius << "\t" << ncables << "\t" << dphi*180/3.14 << G4endl;

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
  //std::cout << "PHG4BbcDetector::ConstructSupport()" << std::endl;

  G4Material *Aluminum = GetDetectorMaterial("G4_Al");

  // BBC Base Plates
  G4double basep_width = 35.56 * cm;
  G4double basep_height = 1.91 * cm;
  G4double basep_len = 46.99 * cm;
  G4Box *bbc_base_plate = new G4Box("bbc_base_plate", basep_width/2, basep_height/2, basep_len/2);
  G4LogicalVolume *bbc_base_plate_lv = new G4LogicalVolume(bbc_base_plate, Aluminum, G4String("Bbc_Base_Plates"));
  GetDisplayAction()->AddVolume(bbc_base_plate_lv, "Bbc_Base_Plates");

  G4VPhysicalVolume *base_plate_vol[2] = {nullptr};  // Mount Plates

  // Place South Base Plates
  base_plate_vol[0] = new G4PVPlacement(nullptr, G4ThreeVector(0, -15.*cm - basep_height/2, (-250. -12.) * cm),
                                    bbc_base_plate_lv, "BBC_BASE_PLATE", logicWorld, false, 0, OverlapCheck());

  // Place North Base Plates
  base_plate_vol[1] = new G4PVPlacement(nullptr, G4ThreeVector(0, -15.*cm - basep_height/2, (250. + 12.0) * cm),
                                    bbc_base_plate_lv, "BBC_BASE_PLATE", logicWorld, false, 1, OverlapCheck());

  // this is more to prevent compiler warnings about unused variables
  if (!base_plate_vol[0] || !base_plate_vol[1])
  {
    std::cout << "Problem placing BBC supports" << std::endl;
  }
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
