#include "PHG4ZDCDetector.h"

#include "PHG4ZDCDefs.h"
#include "PHG4ZDCDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4gdml/PHG4GDMLUtility.hh>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <phool/recoConsts.h>

#include <TSystem.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>            // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>

class G4Material;
class G4VSolid;
class PHCompositeNode;

//_______________________________________________________________________
PHG4ZDCDetector::PHG4ZDCDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam, const int detid)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4ZDCDisplayAction*>(subsys->GetDisplayAction()))
  , m_Params(parameters)
  , m_GdmlConfig(PHG4GDMLUtility::GetOrMakeConfigNode(Node))
  , m_Angle(M_PI_4 * radian)
  , m_TPlate(2.3 * mm)
  , m_HPlate(400.0 * mm)
  , m_WPlate(100.0 * mm)
  , m_TAbsorber(5.0 * mm)
  , m_HAbsorber(150.0 * mm)
  , m_WAbsorber(100.0 * mm)
  , m_DFiber(0.5 * mm)
  , m_HFiber(400.0 * mm)
  , m_WFiber(100.0 * mm)
  , m_GFiber(0.0001 * mm)
  , m_Gap(0.2 * mm)
  , m_TSMD(10.0 * mm)
  , m_HSMD(160.0 * mm)
  , m_WSMD(105.0 * mm)
  , m_RHole(63.9 * mm)
  , m_TWin(4.7 * mm)
  , m_RWin(228.60 * mm)
  , m_PlaceHole(122.56 * mm)
  , m_Pxwin(0.0 * mm)
  , m_Pywin(0.0 * mm)
  , m_Pzwin(915.5 * mm)
  , m_NMod(3)
  , m_NLay(27)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_AbsorberActiveFlag(m_Params->get_int_param("absorberactive"))
  , m_SupportActiveFlag(m_Params->get_int_param("supportactive"))
  , m_Layer(detid)
  , m_SuperDetector("NONE")
{
  assert(m_GdmlConfig);
}

//_______________________________________________________________________
int PHG4ZDCDetector::IsInZDC(G4VPhysicalVolume* volume) const
{
  G4LogicalVolume* mylogvol = volume->GetLogicalVolume();

  if (m_ActiveFlag)
  {
    if (m_ScintiLogicalVolSet.find(mylogvol) != m_ScintiLogicalVolSet.end())
    {
      return 1;
    }
    if (m_FiberLogicalVolSet.find(mylogvol) != m_FiberLogicalVolSet.end())
    {
      return 2;
    }
  }
  if (m_AbsorberActiveFlag)
  {
    if (m_AbsorberLogicalVolSet.find(mylogvol) != m_AbsorberLogicalVolSet.end())
    {
      return -1;
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

//_______________________________________________________________________
void PHG4ZDCDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if (Verbosity() > 0)
  {
    std::cout << "PHG4ZDCDetector: Begin Construction" << std::endl;
  }

  if (m_Layer != PHG4ZDCDefs::NORTH &&
      m_Layer != PHG4ZDCDefs::SOUTH)
  {
    std::cout << "use either PHG4ZDCDefs::NORTH or PHG4ZDCDefs::SOUTH for ZDC Subsystem" << std::endl;
    gSystem->Exit(1);
    return;
  }

  recoConsts* rc = recoConsts::instance();
  G4Material* WorldMaterial = GetDetectorMaterial(rc->get_StringFlag("WorldMaterial"));
  G4Material* Fe = GetDetectorMaterial("G4_Fe");
  G4Material* W = GetDetectorMaterial("G4_W");
  G4Material* PMMA = GetDetectorMaterial("G4_PLEXIGLASS");
  G4Material* Scint = GetDetectorMaterial("G4_POLYSTYRENE");  //place holder

  G4double TGap = m_DFiber + m_Gap;
  G4double Mod_Length = 2 * m_TPlate + m_NLay * (TGap + m_TAbsorber);
  G4double Det_Length = m_NMod * Mod_Length + m_TSMD;
  G4double RTT = 1. / sin(m_Angle);
  G4double First_Pos = -RTT * Det_Length / 2;
  G4double Room = 3.5 * mm;
  G4double Mother_2Z = RTT * Det_Length + 2. * (m_HFiber - m_HAbsorber / 2.) * cos(m_Angle);
  /* Create exit windows */
  G4VSolid* ExitWindow_nocut_solid = new G4Tubs(G4String("ExitWindow_nocut_solid"),
                                                0.0, m_RWin, m_TWin, 0.0, CLHEP::twopi);

  G4VSolid* Hole_solid = new G4Tubs(G4String("Hole_solid"),
                                    0.0, m_RHole, 2 * m_TWin, 0.0, CLHEP::twopi);
  G4VSolid* ExitWindow_1cut_solid = new G4SubtractionSolid("ExitWindow_1cut_solid", ExitWindow_nocut_solid, Hole_solid, nullptr, G4ThreeVector(m_PlaceHole, 0, 0));

  G4VSolid* ExitWindow_2cut_solid = new G4SubtractionSolid("ExitWindow_2cut_solid", ExitWindow_1cut_solid, Hole_solid, nullptr, G4ThreeVector(-m_PlaceHole, 0, 0));

  G4LogicalVolume* ExitWindow_log = new G4LogicalVolume(ExitWindow_2cut_solid, GetDetectorMaterial("G4_STAINLESS-STEEL"), G4String("ExitWindow_log"), nullptr, nullptr, nullptr);

  GetDisplayAction()->AddVolume(ExitWindow_log, "Window");
  m_SupportLogicalVolSet.insert(ExitWindow_log);
  G4RotationMatrix Window_rotm;
  Window_rotm.rotateX(m_Params->get_double_param("rot_x") * deg);
  Window_rotm.rotateY(m_Params->get_double_param("rot_y") * deg);
  Window_rotm.rotateZ(m_Params->get_double_param("rot_z") * deg);

  if (m_Layer == PHG4ZDCDefs::NORTH)
  {
    new G4PVPlacement(G4Transform3D(Window_rotm, G4ThreeVector(m_Pxwin, m_Pywin, m_Params->get_double_param("place_z") * cm - m_Pzwin)),
                      ExitWindow_log, "Window_North", logicWorld, false, PHG4ZDCDefs::NORTH, OverlapCheck());
  }

  else if (m_Layer == PHG4ZDCDefs::SOUTH)
  {
    new G4PVPlacement(G4Transform3D(Window_rotm, G4ThreeVector(m_Pxwin, m_Pywin, -m_Params->get_double_param("place_z") * cm + m_Pzwin)),
                      ExitWindow_log, "Window_South", logicWorld, false, PHG4ZDCDefs::SOUTH, OverlapCheck());
  }
  /* ZDC detector here */
  /* Create the box envelope = 'world volume' for ZDC */

  G4double Mother_X = m_WAbsorber / 2. + Room;
  G4double Mother_Y = (m_HFiber - m_HAbsorber / 2.) * sin(m_Angle) + Room;
  G4double Mother_Z = Mother_2Z / 2. + Room;
  G4VSolid* ZDC_envelope_solid = new G4Box(G4String("ZDC_envelope_solid"),
                                           Mother_X,
                                           Mother_Y,
                                           Mother_Z);

  G4LogicalVolume* ZDC_envelope_log = new G4LogicalVolume(ZDC_envelope_solid, WorldMaterial, G4String("ZDC_envelope_log"), nullptr, nullptr, nullptr);

  /* Define visualization attributes */
  GetDisplayAction()->AddVolume(ZDC_envelope_log, "Envelope");

  /* Define rotation attributes for envelope cone */
  G4RotationMatrix ZDC_rotm;
  ZDC_rotm.rotateX(m_Params->get_double_param("rot_x") * deg);
  ZDC_rotm.rotateY(m_Params->get_double_param("rot_y") * deg);
  ZDC_rotm.rotateZ(m_Params->get_double_param("rot_z") * deg);

  /* Create logical volumes for a plate to contain fibers */
  G4VSolid* fiber_plate_solid = new G4Box(G4String("fiber_plate_solid"),
                                          m_WFiber / 2.,
                                          m_HFiber / 2.,
                                          TGap / 2.);

  G4LogicalVolume* fiber_plate_log = new G4LogicalVolume(fiber_plate_solid, WorldMaterial, G4String("fiber_plate_log"), nullptr, nullptr, nullptr);
  GetDisplayAction()->AddVolume(fiber_plate_log, "fiber_plate_air");
  /*  front and back plate */
  G4VSolid* fb_plate_solid = new G4Box(G4String("fb_plate_solid"),
                                       m_WPlate / 2.,
                                       m_HPlate / 2.,
                                       m_TPlate / 2.);

  G4LogicalVolume* fb_plate_log = new G4LogicalVolume(fb_plate_solid, Fe, G4String("fb_plate_log"), nullptr, nullptr, nullptr);
  m_SupportLogicalVolSet.insert(fb_plate_log);
  GetDisplayAction()->AddVolume(fb_plate_log, "FrontBackPlate");

  /*  absorber */
  G4VSolid* absorber_solid = new G4Box(G4String("absorber_solid"),
                                       m_WAbsorber / 2.,
                                       m_HAbsorber / 2.,
                                       m_TAbsorber / 2.);

  G4LogicalVolume* absorber_log = new G4LogicalVolume(absorber_solid, W, G4String("absorber_log"), nullptr, nullptr, nullptr);
  m_AbsorberLogicalVolSet.insert(absorber_log);
  GetDisplayAction()->AddVolume(absorber_log, "Absorber");

  /* SMD */
  //volume that contains scintillators
  G4VSolid* SMD_solid = new G4Box(G4String("SMD_solid"),
                                  m_WSMD / 2.,
                                  m_HSMD / 2.,
                                  m_TSMD / 2.);

  G4LogicalVolume* SMD_log = new G4LogicalVolume(SMD_solid, WorldMaterial, G4String("SMD_log"), nullptr, nullptr, nullptr);
  GetDisplayAction()->AddVolume(SMD_log, "SMD");
  // small scintillators block
  G4double scintx = 15 * mm;
  G4double scinty = 20 * mm;
  G4VSolid* Scint_solid = new G4Box(G4String("Scint_solid"),
                                    scintx / 2.,
                                    scinty / 2.,
                                    m_TSMD / 2.);

  G4LogicalVolume* Scint_log = new G4LogicalVolume(Scint_solid, Scint, G4String("Scint_log"), nullptr, nullptr, nullptr);
  m_ScintiLogicalVolSet.insert(Scint_log);
  GetDisplayAction()->AddVolume(Scint_log, "Scint_solid");

  //put scintillators in the SMD volume
  double scint_XPos = -m_WSMD / 2.;
  double scint_Xstep = scintx / 2.;
  double scint_Ystep = scinty / 2.;
  int Nx = m_WSMD / scintx;
  int Ny = m_HSMD / scinty;

  for (int i = 0; i < Nx; i++)
  {
    scint_XPos += scint_Xstep;
    double scint_YPos = -m_HSMD / 2.;
    for (int j = 0; j < Ny; j++)
    {
      int copyno = Nx * j + i;
      scint_YPos += scint_Ystep;
      G4RotationMatrix SMDRotation;
      new G4PVPlacement(G4Transform3D(SMDRotation, G4ThreeVector(scint_XPos, scint_YPos, 0.0)),
                        Scint_log,
                        "single_scint",
                        SMD_log,
                        false, copyno, OverlapCheck());

      scint_YPos += scint_Ystep;
    }
    scint_XPos += scint_Xstep;
  }

  /* Create logical volumes for fibers */
  G4VSolid* single_fiber_solid = new G4Tubs(G4String("single_fiber_solid"),
                                            0.0, (m_DFiber / 2.) - m_GFiber, (m_HFiber / 2.), 0.0, CLHEP::twopi);

  G4LogicalVolume* single_fiber_log = new G4LogicalVolume(single_fiber_solid, PMMA, G4String("single_fiber_log"), nullptr, nullptr, nullptr);
  m_FiberLogicalVolSet.insert(single_fiber_log);
  GetDisplayAction()->AddVolume(single_fiber_log, "Fiber");

  /* Rotation Matrix for fibers */
  G4RotationMatrix* FiberRotation = new G4RotationMatrix();
  FiberRotation->rotateX(90. * deg);

  /* Place fibers in the fiber plate */
  G4double fiber_XPos = -m_WFiber / 2.;
  G4double fiber_step = m_DFiber / 2.;
  int Nfiber = m_WFiber / m_DFiber;

  for (int i = 0; i < Nfiber; i++)
  {
    fiber_XPos += fiber_step;
    int copyno = i;

    new G4PVPlacement(FiberRotation, G4ThreeVector(fiber_XPos, 0.0, 0.0),
                      single_fiber_log,
                      G4String("single_fiber_scint"),
                      fiber_plate_log,
                      false, copyno, OverlapCheck());
    fiber_XPos += fiber_step;
  }

  /* Rotation for plates in ZDC */
  G4RotationMatrix* PlateRotation = new G4RotationMatrix();

  PlateRotation->rotateX(-m_Angle);

  /* construct ZDC */
  G4double Plate_Step = m_TPlate * RTT;
  G4double Absorber_Step = m_TAbsorber * RTT;
  G4double Gap_Step = TGap * RTT;
  G4double SMD_Step = m_TSMD * RTT;
  G4double ZPos = First_Pos - Plate_Step / 2.;
  G4double Gap_YPos = (m_HFiber - m_HAbsorber) / 2. * sin(m_Angle);
  G4double Gap_ZPos = (m_HFiber - m_HAbsorber) / 2. * cos(m_Angle);

  G4double Plate_YPos = (m_HPlate - m_HAbsorber) / 2. * sin(m_Angle);
  G4double Plate_ZPos = (m_HPlate - m_HAbsorber) / 2. * cos(m_Angle);

  G4double SMD_YPos = (m_HSMD - m_HAbsorber) / 2. * sin(m_Angle);
  G4double SMD_ZPos = (m_HSMD - m_HAbsorber) / 2. * cos(m_Angle);

  /* start the loop: for every module ---  front plate-absorber-fiber plate-absorber-.....-fiber plate-back plate */
  int copyno_plate = 0;
  for (int i = 0; i < m_NMod; i++)
  {
    //place the SMD in between the 1st and 2nd module
    if (i == 1)
    {
      ZPos += (SMD_Step / 2.);
      new G4PVPlacement(PlateRotation, G4ThreeVector(0.0, SMD_YPos, SMD_ZPos + ZPos),
                        SMD_log,
                        G4String("SMD"),
                        ZDC_envelope_log,
                        false, 0, OverlapCheck());  //using copy number for now, need to find a better way
      ZPos += (SMD_Step / 2.);
    }
    /* place the front plate */
    ZPos += (Plate_Step / 2.);
    new G4PVPlacement(PlateRotation, G4ThreeVector(0.0, Plate_YPos, Plate_ZPos + ZPos),
                      fb_plate_log,
                      G4String("front_plate"),
                      ZDC_envelope_log,
                      false, copyno_plate, OverlapCheck());
    ZPos += (Plate_Step / 2.);
    copyno_plate++;
    for (int j = 0; j < m_NLay; j++)
    {
      /* place the Absorber */
      ZPos += (Absorber_Step / 2.);
      new G4PVPlacement(PlateRotation, G4ThreeVector(0.0, 0.0, ZPos),
                        absorber_log,
                        G4String("single_absorber"),
                        ZDC_envelope_log,
                        false, i * 100 + j, OverlapCheck());
      ZPos += (Absorber_Step / 2.);

      /* place the fiber plate */
      ZPos += (Gap_Step / 2.);
      int copyno = 27 * i + j;
      std::string name_fiber_plate = "Fiber_Plate_" + std::to_string(copyno);
      new G4PVPlacement(PlateRotation, G4ThreeVector(0.0, Gap_YPos, Gap_ZPos + ZPos),
                        fiber_plate_log,
                        name_fiber_plate,
                        ZDC_envelope_log,
                        false, copyno, OverlapCheck());
      ZPos += (Gap_Step / 2.);
    }
    /* place the back plate */
    ZPos += (Plate_Step / 2.);
    new G4PVPlacement(PlateRotation, G4ThreeVector(0.0, Plate_YPos, Plate_ZPos + ZPos),
                      fb_plate_log,
                      G4String("back_plate"),
                      ZDC_envelope_log,
                      false, copyno_plate, OverlapCheck());
    copyno_plate++;
    ZPos += (Plate_Step / 2.);
  }

  /* Place envelope cone in simulation */

  if (m_Layer == PHG4ZDCDefs::NORTH)
  {
    new G4PVPlacement(G4Transform3D(ZDC_rotm, G4ThreeVector(m_Params->get_double_param("place_x") * cm, m_Params->get_double_param("place_y") * cm, m_Params->get_double_param("place_z") * cm)),
                      ZDC_envelope_log, "ZDC_Envelope_North", logicWorld, false, PHG4ZDCDefs::NORTH, OverlapCheck());
  }
  else if (m_Layer == PHG4ZDCDefs::SOUTH)
  {
    ZDC_rotm.rotateY(180 * deg);
    new G4PVPlacement(G4Transform3D(ZDC_rotm, G4ThreeVector(m_Params->get_double_param("place_x") * cm, m_Params->get_double_param("place_y") * cm, -m_Params->get_double_param("place_z") * cm)),
                      ZDC_envelope_log, "ZDC_Envelope_South", logicWorld, false, PHG4ZDCDefs::SOUTH, OverlapCheck());
  }
  return;
}
