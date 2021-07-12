#include "PHG4ZDCDetector.h"

#include "PHG4ZDCDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4gdml/PHG4GDMLConfig.hh>
#include <g4gdml/PHG4GDMLUtility.hh>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <phool/recoConsts.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>            // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume

#include <TSystem.h>

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>  // for pair, make_pair

class G4VSolid;
class PHCompositeNode;

using namespace std;

//_______________________________________________________________________
PHG4ZDCDetector::PHG4ZDCDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4ZDCDisplayAction*>(subsys->GetDisplayAction()))
  , m_Params(parameters)
  , m_GdmlConfig(PHG4GDMLUtility::GetOrMakeConfigNode(Node))
  , m_Angle( M_PI_4 * radian)
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
  , m_XRot(0.0)
  , m_YRot(0.0)
  , m_ZRot(0.0)
  , m_PlaceX(0.0 * mm)
  , m_PlaceY(0.0 * mm)
  , m_PlaceZ(18430.0 * mm)
  , m_NMod(3)
  , m_NLay(27)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_AbsorberActiveFlag(m_Params->get_int_param("absorberactive"))
  , m_Layer(0)
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
    }
  if (m_AbsorberActiveFlag)
    {
      if (m_AbsorberLogicalVolSet.find(mylogvol) != m_AbsorberLogicalVolSet.end())
	{
	  return -1;
	}
    }
  return 0;
}


//_______________________________________________________________________
void PHG4ZDCDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if ( Verbosity() > 0 )
  {
    cout << "PHG4ZDCDetector: Begin Construction" << endl;
  }
  
  recoConsts* rc = recoConsts::instance();
  G4Material* WorldMaterial = G4Material::GetMaterial(rc->get_StringFlag("WorldMaterial"));
  G4Material* Fe = G4Material::GetMaterial("G4_Fe");
  G4Material* W = G4Material::GetMaterial("G4_W");
  G4Material* PMMA = G4Material::GetMaterial("G4_PLEXIGLASS");
    
  G4double TGap = m_DFiber+m_Gap;
  G4double Mod_Length = 2 * m_TPlate + m_NLay * (TGap + m_TAbsorber);
  G4double Det_Length = m_NMod * Mod_Length;
  G4double RTT = 1./sin(m_Angle);
  G4double First_Pos = -RTT * Det_Length / 2;
  G4double Room = 10 * mm;
  G4double Mother_2Z = RTT * Det_Length + 2. * (m_HFiber - m_HAbsorber / 2.) * cos(m_Angle);
  
  /* Create the box envelope = 'world volume' for ZDC */
  
  G4double Mother_X = m_WAbsorber / 2. + Room;
  G4double Mother_Y = (m_HFiber - m_HAbsorber / 2.) * sin(m_Angle) + Room;
  G4double Mother_Z = Mother_2Z / 2. + Room;
  G4VSolid* ZDC_envelope_solid = new G4Box(G4String("ZDC_envelope_solid"),
                                             Mother_X,
                                             Mother_Y,
                                             Mother_Z);

  G4LogicalVolume* ZDC_envelope_log = new G4LogicalVolume(ZDC_envelope_solid, WorldMaterial, G4String("ZDC_envelope_log"), 0, 0, 0);

  /* Define visualization attributes */
  GetDisplayAction()->AddVolume(ZDC_envelope_log, "Envelope");

  /* Define rotation attributes for envelope cone */
  G4RotationMatrix ZDC_rotm;
  ZDC_rotm.rotateX(m_XRot);
  ZDC_rotm.rotateY(m_YRot);
  ZDC_rotm.rotateZ(m_ZRot);

  /* Place envelope cone in simulation */
  string name_envelope = "ZDC_phy_envelope_0";
  
  new G4PVPlacement(G4Transform3D(ZDC_rotm, G4ThreeVector(m_PlaceX, m_PlaceY, m_PlaceZ)),
                    ZDC_envelope_log, name_envelope, logicWorld, 0, 0, OverlapCheck());

  /* Create logical volumes for a plate to contain fibers */
  G4VSolid* fiber_plate_solid = new G4Box(G4String("fiber_plate_solid"),
					  m_WFiber / 2.,
					  m_HFiber / 2.,
					  TGap / 2.);
  
  G4LogicalVolume* fiber_plate_log = new G4LogicalVolume(fiber_plate_solid, WorldMaterial, G4String("fiber_plate_log"), 0, 0, 0);					  
  m_AbsorberLogicalVolSet.insert(fiber_plate_log);

  /*  front and back plate */
  G4VSolid* fb_plate_solid = new G4Box(G4String("fb_plate_solid"),
				       m_WPlate / 2.,
				       m_HPlate / 2.,
				       m_TPlate / 2.);

  G4LogicalVolume* fb_plate_log = new G4LogicalVolume(fb_plate_solid, Fe, G4String("fb_plate_log"), 0, 0, 0);
  m_AbsorberLogicalVolSet.insert(fb_plate_log);
  GetDisplayAction()->AddVolume(fb_plate_log, "FrontBackPlate");

  /*  absorber */
  G4VSolid* absorber_solid = new G4Box(G4String("absorber_solid"),
				       m_WAbsorber / 2.,
				       m_HAbsorber / 2.,
				       m_TAbsorber / 2.);
  
  G4LogicalVolume* absorber_log = new G4LogicalVolume(absorber_solid, W, G4String("absorber_log"), 0, 0, 0);
  m_AbsorberLogicalVolSet.insert(absorber_log);
  GetDisplayAction()->AddVolume(absorber_log, "Absorber");

  /* Create logical volumes for fibers */
  G4VSolid* single_fiber_solid = new G4Tubs(G4String("single_fiber_solid"),
					    0.0, (m_DFiber / 2.) - m_GFiber, (m_HFiber / 2.), 0.0, CLHEP::twopi);
  
  G4LogicalVolume* single_fiber_log = new G4LogicalVolume(single_fiber_solid, PMMA, G4String("single_fiber_log"), 0, 0, 0);
  m_ScintiLogicalVolSet.insert(single_fiber_log);
  GetDisplayAction()->AddVolume(single_fiber_log, "Fiber");

  
  /* Rotation Matrix for fibers */
  G4RotationMatrix* FiberRotation = new G4RotationMatrix();
  FiberRotation->rotateX(90.*deg);

  /* Place fibers in the fiber plate */
  G4double fiber_XPos = -m_WFiber / 2.;
  G4double fiber_step = m_DFiber / 2.;
  int Nfiber = m_WFiber / m_DFiber;
  
  for(int i=0; i<Nfiber; i++)
    {
      fiber_XPos+=fiber_step;
    
      new G4PVPlacement(FiberRotation, G4ThreeVector(fiber_XPos, 0.0, 0.0),
                        single_fiber_log,
			G4String("single_fiber_scint"),
                        fiber_plate_log,
                        0, 0, OverlapCheck());
      fiber_XPos+=fiber_step;
    }
  GetDisplayAction()->AddVolume(fiber_plate_log, "FiberPlate");

  /* Rotation for plates in ZDC */
  G4RotationMatrix* PlateRotation= new G4RotationMatrix();;
  PlateRotation->rotateX( - m_Angle);
  
   /* construct ZDC */
  G4double Plate_Step = m_TPlate * RTT;
  G4double Absorber_Step = m_TAbsorber * RTT;
  G4double Gap_Step = TGap * RTT;
  G4double ZPos = First_Pos - Plate_Step / 2.;
  G4double Gap_YPos = (m_HFiber - m_HAbsorber) / 2. * sin(m_Angle);
  G4double Gap_ZPos = (m_HFiber - m_HAbsorber) / 2. * cos(m_Angle);

  G4double Plate_YPos = (m_HPlate - m_HAbsorber) / 2. * sin(m_Angle);
  G4double Plate_ZPos = (m_HPlate - m_HAbsorber) / 2. * cos(m_Angle);
  
  /* start the loop: for every module ---  front plate-absorber-fiber plate-absorber-.....-fiber plate-back plate */
  for(int i=0; i<m_NMod; i++)
    {
      
      /* place the front plate */
      ZPos += (Plate_Step / 2.);
      new G4PVPlacement(PlateRotation, G4ThreeVector(0.0, Plate_YPos, Plate_ZPos + ZPos),
                        fb_plate_log,
                        G4String("front_plate"),
                        ZDC_envelope_log,
                        0, 0, OverlapCheck());
      ZPos += (Plate_Step / 2.);
      for(int j=0; j<m_NLay; j++)
	{

	  /* place the Absorber */
	  ZPos += (Absorber_Step / 2.);
	  new G4PVPlacement(PlateRotation, G4ThreeVector(0.0, 0.0, ZPos),
			    absorber_log,
			    G4String("single_absorber"),
			    ZDC_envelope_log,
			    0, 0, OverlapCheck());
	  ZPos += (Absorber_Step / 2.);
	  
	  /* place the fiber plate */
	  ZPos += (Gap_Step / 2.);
	  ostringstream name_fiber_plate;
	  int copyno = 27 * i + j;
	  name_fiber_plate <<  "Fiber_Plate_" << copyno;
	  new G4PVPlacement(PlateRotation, G4ThreeVector(0.0, Gap_YPos, Gap_ZPos + ZPos),
			    fiber_plate_log,
			    name_fiber_plate.str().c_str(),
			    ZDC_envelope_log,
			    0, copyno, OverlapCheck());
	  ZPos += (Gap_Step / 2.);
	}
      /* place the back plate */
      ZPos += (Plate_Step / 2.);
      new G4PVPlacement(PlateRotation, G4ThreeVector(0.0, Plate_YPos, Plate_ZPos + ZPos),
                        fb_plate_log,
                        G4String("back_plate"),
                        ZDC_envelope_log,
                        0, 0, OverlapCheck());
      ZPos += (Plate_Step / 2.);
    }
  /* place the other ZDC */
  ZDC_rotm.rotateY(180 * deg);
  name_envelope = "ZDC_phy_envelope_1";
  new G4PVPlacement(G4Transform3D(ZDC_rotm, G4ThreeVector(m_PlaceX, m_PlaceY, -m_PlaceZ)),
  ZDC_envelope_log, name_envelope, logicWorld, 0, 1, OverlapCheck());

  return;
 }
