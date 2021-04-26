#include "PHG4CylinderDetector.h"
#include "PHG4CylinderDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <phool/phool.h>

#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4UserLimits.hh>

#include <TSystem.h>

#include <cmath>
#include <iostream>  // for operator<<, endl, basic_ost...
#include <sstream>

class G4VSolid;
class PHCompositeNode;

using namespace std;

//_______________________________________________________________
PHG4CylinderDetector::PHG4CylinderDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int lyr)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(parameters)
  , m_CylinderPhysicalVolume(nullptr)
  , m_DisplayAction(dynamic_cast<PHG4CylinderDisplayAction *>(subsys->GetDisplayAction()))
  , m_Layer(lyr)
{
}

//_______________________________________________________________
bool PHG4CylinderDetector::IsInCylinder(const G4VPhysicalVolume *volume) const
{
  if (volume == m_CylinderPhysicalVolume)
  {
    return true;
  }
  return false;
}

//_______________________________________________________________
void PHG4CylinderDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  G4Material *TrackerMaterial = G4Material::GetMaterial(m_Params->get_string_param("material"));

  if (!TrackerMaterial)
  {
    cout << "Error: Can not set material " << m_Params->get_string_param("material") << endl;
    gSystem->Exit(1);
  }

  // determine length of cylinder using PHENIX's rapidity coverage if flag is true
  double radius = m_Params->get_double_param("radius") * cm;
  double thickness = m_Params->get_double_param("thickness") * cm;
  double length = m_Params->get_double_param("length") * cm;
  if (!isfinite(radius) || !isfinite(thickness) || !isfinite(length))
  {
    cout << PHWHERE << ": Bad Parameters for " << GetName() << endl;
    cout << "Radius: " << radius << endl;
    cout << "Thickness: " << thickness << endl;
    cout << "Length: " << length << endl;
    gSystem->Exit(1);
  }
  G4VSolid *cylinder_solid = new G4Tubs(G4String(GetName()),
                                        radius,
                                        radius + thickness,
                                        length / 2., 0, twopi);
  double steplimits = m_Params->get_double_param("steplimits") * cm;
  G4UserLimits *g4userlimits = nullptr;
  if (isfinite(steplimits))
  {
    g4userlimits = new G4UserLimits(steplimits);
  }

  G4LogicalVolume *cylinder_logic = new G4LogicalVolume(cylinder_solid,
                                                        TrackerMaterial,
                                                        G4String(GetName()),
                                                        nullptr, nullptr, g4userlimits);
  PHG4Subsystem *mysys = GetMySubsystem();
  mysys->SetLogicalVolume(cylinder_logic);

  G4RotationMatrix *rotm = new G4RotationMatrix();
  rotm->rotateY(m_Params->get_double_param("rot_y") * deg);

  m_CylinderPhysicalVolume = new G4PVPlacement(rotm,
					       G4ThreeVector(m_Params->get_double_param("place_x") * cm,
							     m_Params->get_double_param("place_y") * cm,
							     m_Params->get_double_param("place_z") * cm),
                                               cylinder_logic,
                                               G4String(GetName()),
                                               logicWorld, 0, false, OverlapCheck());
  m_DisplayAction->SetMyVolume(cylinder_logic);
}
