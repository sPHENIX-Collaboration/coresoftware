#include "PHG4CylinderDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <sstream>

using namespace std;

//_______________________________________________________________
PHG4CylinderDetector::PHG4CylinderDetector(PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int lyr)
  : PHG4Detector(Node, dnam)
  , params(parameters)
  , cylinder_physi(nullptr)
  , layer(lyr)
{
}

//_______________________________________________________________
bool PHG4CylinderDetector::IsInCylinder(const G4VPhysicalVolume *volume) const
{
  if (volume == cylinder_physi)
  {
    return true;
  }
  return false;
}

//_______________________________________________________________
void PHG4CylinderDetector::Construct(G4LogicalVolume *logicWorld)
{
  G4Material *TrackerMaterial = G4Material::GetMaterial(params->get_string_param("material"));

  if (!TrackerMaterial)
  {
    std::cout << "Error: Can not set material" << std::endl;
    exit(-1);
  }

  G4VisAttributes *siliconVis = new G4VisAttributes();
  if (params->get_int_param("blackhole"))
  {
    PHG4Utils::SetColour(siliconVis, "BlackHole");
    siliconVis->SetVisibility(false);
    siliconVis->SetForceSolid(false);
  }
  else
  {
    PHG4Utils::SetColour(siliconVis, params->get_string_param("material"));
    siliconVis->SetVisibility(true);
    siliconVis->SetForceSolid(true);
  }

  // determine length of cylinder using PHENIX's rapidity coverage if flag is true
  double radius = params->get_double_param("radius") * cm;
  double thickness = params->get_double_param("thickness") * cm;
  G4VSolid *cylinder_solid = new G4Tubs(G4String(GetName()),
                                        radius,
                                        radius + thickness,
                                        params->get_double_param("length") * cm / 2., 0, twopi);
  double steplimits = params->get_double_param("steplimits") * cm;
  G4UserLimits *g4userlimits = nullptr;
  if (isfinite(steplimits))
  {
    g4userlimits = new G4UserLimits(steplimits);
  }

  G4LogicalVolume *cylinder_logic = new G4LogicalVolume(cylinder_solid,
                                                        TrackerMaterial,
                                                        G4String(GetName()),
                                                        nullptr, nullptr, g4userlimits);
  cylinder_logic->SetVisAttributes(siliconVis);
  cylinder_physi = new G4PVPlacement(0, G4ThreeVector(params->get_double_param("place_x") * cm,
                                                      params->get_double_param("place_y") * cm,
                                                      params->get_double_param("place_z") * cm),
                                     cylinder_logic,
                                     G4String(GetName()),
                                     logicWorld, 0, false, OverlapCheck());
}
