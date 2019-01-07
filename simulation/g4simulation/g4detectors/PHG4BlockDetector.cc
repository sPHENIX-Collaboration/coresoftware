#include "PHG4BlockDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4VisAttributes.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________
PHG4BlockDetector::PHG4BlockDetector(PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int lyr)
  : PHG4Detector(Node, dnam)
  , m_Params(parameters)
  , m_BlockPhysi(nullptr)
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
void PHG4BlockDetector::Construct(G4LogicalVolume *logicWorld)
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
  G4VisAttributes matVis;
  if (m_Params->get_int_param("blackhole"))
  {
    PHG4Utils::SetColour(&matVis, "BlackHole");
    matVis.SetVisibility(false);
    matVis.SetForceSolid(false);
  }
  else
  {
    PHG4Utils::SetColour(&matVis, m_Params->get_string_param("material"));
    matVis.SetVisibility(true);
    matVis.SetForceSolid(true);
  }
  block_logic->SetVisAttributes(matVis);

  G4RotationMatrix *rotm = new G4RotationMatrix();
  rotm->rotateZ(m_Params->get_double_param("rot_z") * deg);
  m_BlockPhysi = new G4PVPlacement(rotm, G4ThreeVector(m_Params->get_double_param("place_x") * cm, m_Params->get_double_param("place_y") * cm, m_Params->get_double_param("place_z") * cm),
                                   block_logic,
                                   G4String(GetName()),
                                   logicWorld, 0, false, OverlapCheck());
}
