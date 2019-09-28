#include "PHG4BlockDetector.h"

#include "PHG4BlockDisplayAction.h"
#include "PHG4BlockSubsystem.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>         // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>    // for PHG4DisplayAction

#include <Geant4/G4Box.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>    // for G4RotationMatrix
#include <Geant4/G4String.hh>            // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>       // for G4ThreeVector
#include <Geant4/G4UserLimits.hh>

#include <CLHEP/Units/SystemOfUnits.h>   // for cm, deg

#include <cmath>                         // for isfinite
#include <cstdlib>                      // for exit
#include <iostream>                      // for operator<<, endl, basic_ostream
#include <sstream>

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
  G4Material *TrackerMaterial = G4Material::GetMaterial(m_Params->get_string_param("material"));

  if (!TrackerMaterial)
  {
    std::cout << "Error: Can not set material" << std::endl;
    exit(-1);
  }

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

  G4RotationMatrix *rotm = new G4RotationMatrix();
  rotm->rotateZ(m_Params->get_double_param("rot_z") * deg);
  m_BlockPhysi = new G4PVPlacement(rotm, G4ThreeVector(m_Params->get_double_param("place_x") * cm, m_Params->get_double_param("place_y") * cm, m_Params->get_double_param("place_z") * cm),
                                   block_logic,
                                   G4String(GetName()),
                                   logicWorld, 0, false, OverlapCheck());
  m_DisplayAction->SetMyVolume(block_logic);
}
