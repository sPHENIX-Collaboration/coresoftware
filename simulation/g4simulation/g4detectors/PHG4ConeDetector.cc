#include "PHG4ConeDetector.h"
#include "PHG4ConeDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>         // for PHG4Detector
#include <g4main/PHG4Subsystem.h>
#include <g4main/PHG4Utils.h>

#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>            // for G4String
#include <Geant4/G4SystemOfUnits.hh>   // for cm
#include <Geant4/G4ThreeVector.hh>     // for G4ThreeVector
#include <Geant4/G4VisAttributes.hh>

#include <cmath>     // for M_PI
#include <cstdlib>   // for exit
#include <iostream>  // for operator<<, endl, basic_ostream
#include <sstream>

using namespace std;

//_______________________________________________________________
//note this inactive thickness is ~1.5% of a radiation length
PHG4ConeDetector::PHG4ConeDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int lyr)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(parameters)
  , m_DisplayAction(dynamic_cast<PHG4ConeDisplayAction *>(subsys->GetDisplayAction()))
  , layer(lyr)
{
}

//_______________________________________________________________
//_______________________________________________________________
bool PHG4ConeDetector::IsInConeActive(G4VPhysicalVolume *volume)
{
  if (volume == m_ConePhysVol)
  {
    return true;
  }
  return false;
}

//_______________________________________________________________
void PHG4ConeDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  G4Material *TrackerMaterial = G4Material::GetMaterial(m_Params->get_string_param("material"));

  if (!TrackerMaterial)
  {
    std::cout << "Error: Can not set material" << std::endl;
    exit(-1);
  }

  G4VSolid *cone_solid = new G4Cons((GetName()+"_SOLID"),
                           m_Params->get_double_param("rmin1")*cm,
			   m_Params->get_double_param("rmax1")*cm,
                           m_Params->get_double_param("rmin2")*cm,
			   m_Params->get_double_param("rmax2")*cm,
			   m_Params->get_double_param("length")*cm,
			   m_Params->get_double_param("sphi")*deg,
			   m_Params->get_double_param("dphi")*deg);

  G4LogicalVolume *cone_logic = new G4LogicalVolume(cone_solid,
                                    TrackerMaterial,
                                    G4String(GetName()+"_LOGIC"),
                                    0, 0, 0);
  PHG4Subsystem *mysys = GetMySubsystem();
  mysys->SetLogicalVolume(cone_logic);

  G4RotationMatrix *rotm = new G4RotationMatrix();
  rotm->rotateZ( m_Params->get_double_param("rot_z")*deg);
  m_ConePhysVol = new G4PVPlacement(rotm, 
				    G4ThreeVector(m_Params->get_double_param("place_x") * cm,
						  m_Params->get_double_param("place_y") * cm,
						  m_Params->get_double_param("place_z") * cm),
				    cone_logic,
				    G4String(GetName()),
				    logicWorld, 0, false, OverlapCheck());
  m_DisplayAction->SetMyVolume(cone_logic);
}
