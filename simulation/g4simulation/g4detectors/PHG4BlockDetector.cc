#include "PHG4BlockDetector.h"

#include "PHG4BlockDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4UserLimits.hh>

#include <CLHEP/Units/SystemOfUnits.h>  // for cm, deg

#include <cmath>     // for isfinite
#include <iostream>  // for operator<<, endl, basic_ostream
#include <sstream>

class G4Material;
class G4VSolid;
class PHCompositeNode;

using namespace std;

//_______________________________________________________________
PHG4BlockDetector::PHG4BlockDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int lyr)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(parameters)
  , m_BlockPhysi(nullptr)
  , m_DisplayAction(dynamic_cast<PHG4BlockDisplayAction *>(subsys->GetDisplayAction()))
  , m_Layer(lyr)
{
}

//_______________________________________________________________
bool PHG4BlockDetector::IsInBlock(G4VPhysicalVolume *volume) const
{
  if (volume == m_BlockPhysi)
  {
    return true;
  }
  return false;
}

//_______________________________________________________________
void PHG4BlockDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  G4Material *TrackerMaterial = GetDetectorMaterial(m_Params->get_string_param("material"));

  G4VSolid *block_solid = new G4Box(G4String(GetName()),
                                    m_Params->get_double_param("size_x") / 2. * cm,
                                    m_Params->get_double_param("size_y") / 2. * cm,
                                    m_Params->get_double_param("size_z") / 2. * cm);

  double steplimits = m_Params->get_double_param("steplimits") * cm;
  G4UserLimits *g4userlimits = nullptr;
  if (isfinite(steplimits))
  {
    g4userlimits = new G4UserLimits(steplimits);
  }

  G4LogicalVolume *block_logic = new G4LogicalVolume(block_solid,
                                                     TrackerMaterial,
                                                     G4String(GetName()),
                                                     nullptr, nullptr, g4userlimits);

  PHG4Subsystem *mysys = GetMySubsystem();
  mysys->SetLogicalVolume(block_logic);

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
    cout << __PRETTY_FUNCTION__ << ": Warning : " << GetName() << " is configured with more than one of the x-y-z rotations of "
         << "(" << m_Params->get_double_param("rot_x") << ", "
         << m_Params->get_double_param("rot_x") << ", "
         << m_Params->get_double_param("rot_x") << ") degrees. "
         << "The rotation is instruction is ambiguous and they are performed in the order of X->Y->Z rotations with result rotation matrix of:";
    rotm->print(cout);
  }

  m_BlockPhysi = new G4PVPlacement(rotm, G4ThreeVector(m_Params->get_double_param("place_x") * cm, m_Params->get_double_param("place_y") * cm, m_Params->get_double_param("place_z") * cm),
                                   block_logic,
                                   G4String(GetName()),
                                   logicWorld, false, false, OverlapCheck());
  m_DisplayAction->SetMyVolume(block_logic);
}
