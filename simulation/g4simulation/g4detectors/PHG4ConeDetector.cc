#include "PHG4ConeDetector.h"
#include "PHG4ConeDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SystemOfUnits.hh>   // for cm
#include <Geant4/G4ThreeVector.hh>     // for G4ThreeVector

#include <iostream>  // for operator<<, endl, basic_ostream
#include <sstream>

class G4Material;
class G4VSolid;

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
  G4Material *TrackerMaterial = GetDetectorMaterial(m_Params->get_string_param("material"));

  G4VSolid *cone_solid = new G4Cons((GetName() + "_SOLID"),
                                    m_Params->get_double_param("rmin1") * cm,
                                    m_Params->get_double_param("rmax1") * cm,
                                    m_Params->get_double_param("rmin2") * cm,
                                    m_Params->get_double_param("rmax2") * cm,
                                    m_Params->get_double_param("length") * cm,
                                    m_Params->get_double_param("sphi") * deg,
                                    m_Params->get_double_param("dphi") * deg);

  G4LogicalVolume *cone_logic = new G4LogicalVolume(cone_solid,
                                                    TrackerMaterial,
                                                    GetName() + "_LOGIC",
                                                    nullptr, nullptr, nullptr);
  PHG4Subsystem *mysys = GetMySubsystem();
  mysys->SetLogicalVolume(cone_logic);

  G4RotationMatrix *rotm = new G4RotationMatrix();

  int nRotation(0);
  if (m_Params->get_double_param("rot_x") != 0)
  {
    ++nRotation;
    rotm->rotateX(m_Params->get_double_param("rot_x") * deg);
  }
  if (m_Params->get_double_param("rot_y") != 0)
  {
    ++nRotation;
    rotm->rotateY(m_Params->get_double_param("rot_y") * deg);
  }
  if (m_Params->get_double_param("rot_z") != 0)
  {
    ++nRotation;
    rotm->rotateZ(m_Params->get_double_param("rot_z") * deg);
  }

  if (nRotation >= 2)
  {
    std::cout << __PRETTY_FUNCTION__ << ": Warning : " << GetName() << " is configured with more than one of the x-y-z rotations of "
              << "(" << m_Params->get_double_param("rot_x") << ", "
              << m_Params->get_double_param("rot_x") << ", "
              << m_Params->get_double_param("rot_x") << ") degrees. "
              << "The rotation is instruction is ambiguous and they are performed in the order of X->Y->Z rotations with result rotation matrix of:";
    rotm->print(std::cout);
  }

  m_ConePhysVol = new G4PVPlacement(rotm,
                                    G4ThreeVector(m_Params->get_double_param("place_x") * cm,
                                                  m_Params->get_double_param("place_y") * cm,
                                                  m_Params->get_double_param("place_z") * cm),
                                    cone_logic,
                                    GetName(),
                                    logicWorld, false, false, OverlapCheck());
  m_DisplayAction->SetMyVolume(cone_logic);
}
