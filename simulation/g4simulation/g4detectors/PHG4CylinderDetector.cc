#include "PHG4CylinderDetector.h"
#include "PHG4CylinderDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <phool/phool.h>

#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4UserLimits.hh>

#include <TSystem.h>

#include <cmath>
#include <iostream>  // for operator<<, endl, basic_ost...
#include <sstream>

class G4Material;
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
  G4Material *TrackerMaterial = GetDetectorMaterial(m_Params->get_string_param("material"));

  // determine length of cylinder using PHENIX's rapidity coverage if flag is true
  double radius = m_Params->get_double_param("radius") * cm;
  double thickness = m_Params->get_double_param("thickness") * cm;
  double length = m_Params->get_double_param("length") * cm;

  double start_phi_rad = m_Params->get_double_param("start_phi_rad") * rad;
  double delta_phi_rad = m_Params->get_double_param("delta_phi_rad") * rad;

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
                                        length / 2., start_phi_rad, delta_phi_rad);
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

  m_CylinderPhysicalVolume = new G4PVPlacement(rotm,
                                               G4ThreeVector(m_Params->get_double_param("place_x") * cm,
                                                             m_Params->get_double_param("place_y") * cm,
                                                             m_Params->get_double_param("place_z") * cm),
                                               cylinder_logic,
                                               G4String(GetName()),
                                               logicWorld, false, false, OverlapCheck());
  m_DisplayAction->SetMyVolume(cylinder_logic);
}
