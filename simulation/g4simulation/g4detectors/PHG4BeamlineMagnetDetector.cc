#include "PHG4BeamlineMagnetDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>  // for PHG4Detector
#include <g4main/PHG4Utils.h>

#include <phool/phool.h>

#include <Geant4/G4ChordFinder.hh>
#include <Geant4/G4ClassicalRK4.hh>
#include <Geant4/G4FieldManager.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Mag_UsualEqRhs.hh>
#include <Geant4/G4MagneticField.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4QuadrupoleMagField.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>  // for G4double, G4bool
#include <Geant4/G4UniformMagField.hh>
#include <Geant4/G4VisAttributes.hh>

#include <CLHEP/Units/SystemOfUnits.h>  // for cm, deg, tesla, twopi, meter

#include <cstdlib>   // for exit
#include <iostream>  // for operator<<, basic_ostream

class G4Material;
class G4VSolid;
class PHCompositeNode;
class PHG4Subsystem;

using namespace std;

//_______________________________________________________________
PHG4BeamlineMagnetDetector::PHG4BeamlineMagnetDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int lyr)
  : PHG4Detector(subsys, Node, dnam)
  , params(parameters)
  , magnet_physi(nullptr)
  , cylinder_physi(nullptr)
  , layer(lyr)
{
}

//_______________________________________________________________
bool PHG4BeamlineMagnetDetector::IsInBeamlineMagnet(const G4VPhysicalVolume *volume) const
{
  if (volume == magnet_physi)
  {
    return true;
  }
  return false;
}

//_______________________________________________________________
void PHG4BeamlineMagnetDetector::ConstructMe(G4LogicalVolume *logicMother)
{
  G4Material *TrackerMaterial = GetDetectorMaterial(params->get_string_param("material"));

  G4VisAttributes *fieldVis = new G4VisAttributes();
  PHG4Utils::SetColour(fieldVis, "BlackHole");
  fieldVis->SetVisibility(false);
  fieldVis->SetForceSolid(false);

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

  /* Define origin vector (center of magnet) */
  G4ThreeVector origin(params->get_double_param("place_x") * cm,
                       params->get_double_param("place_y") * cm,
                       params->get_double_param("place_z") * cm);

  /* Define magnet rotation matrix */
  G4RotationMatrix *rotm = new G4RotationMatrix();
  rotm->rotateX(params->get_double_param("rot_x") * deg);
  rotm->rotateY(params->get_double_param("rot_y") * deg);
  rotm->rotateZ(params->get_double_param("rot_z") * deg);

  /* Creating a magnetic field */
  G4MagneticField *magField = nullptr;

  string magnettype = params->get_string_param("magtype");
  if (magnettype == "dipole")
  {
    G4double fieldValue = params->get_double_param("field_y") * tesla;
    magField = new G4UniformMagField(G4ThreeVector(0., fieldValue, 0.));

    if (Verbosity() > 0)
      cout << "Creating DIPOLE with field " << fieldValue << " and name " << GetName() << endl;
  }
  else if (magnettype == "quadrupole")
  {
    G4double fieldGradient = params->get_double_param("fieldgradient") * tesla / meter;

    /* G4MagneticField::GetFieldValue( pos*, B* ) uses GLOBAL coordinates, not local.
       * Therefore, place magnetic field center at the correct location and angle for the
       * magnet AND do the same transformations for the logical volume (see below). */
    magField = new G4QuadrupoleMagField(fieldGradient, origin, rotm);
    //      magField = new PHG4QuadrupoleMagField ( fieldGradient, origin, rotm );

    if (Verbosity() > 0)
    {
      cout << "Creating QUADRUPOLE with gradient " << fieldGradient << " and name " << GetName() << endl;
      cout << "at x, y, z = " << origin.x() << " , " << origin.y() << " , " << origin.z() << endl;
      cout << "with rotation around x, y, z axis  by: " << rotm->phiX() << ", " << rotm->phiY() << ", " << rotm->phiZ() << endl;
    }
  }

  if (!magField)
  {
    cout << PHWHERE << " Aborting, no magnetic field specified for " << GetName() << endl;
    exit(1);
  }

  /* Set up Geant4 field manager */
  G4Mag_UsualEqRhs *localEquation = new G4Mag_UsualEqRhs(magField);
  G4ClassicalRK4 *localStepper = new G4ClassicalRK4(localEquation);
  G4double minStep = 0.25 * mm;  // minimal step, 1 mm is default
  G4ChordFinder *localChordFinder = new G4ChordFinder(magField, minStep, localStepper);

  G4FieldManager *fieldMgr = new G4FieldManager();
  fieldMgr->SetDetectorField(magField);
  fieldMgr->SetChordFinder(localChordFinder);

  /* Add volume with magnetic field */
  double radius = params->get_double_param("radius") * cm;
  double thickness = params->get_double_param("thickness") * cm;

  G4VSolid *magnet_solid = new G4Tubs(GetName(),
                                      0,
                                      radius + thickness,
                                      params->get_double_param("length") * cm / 2., 0, twopi);

  G4LogicalVolume *magnet_logic = new G4LogicalVolume(magnet_solid,
                                                      GetDetectorMaterial("G4_Galactic"),
                                                      GetName(),
                                                      nullptr, nullptr, nullptr);
  magnet_logic->SetVisAttributes(fieldVis);

  /* Set field manager for logical volume */
  G4bool allLocal = true;
  magnet_logic->SetFieldManager(fieldMgr, allLocal);

  /* create magnet physical volume */
  magnet_physi = new G4PVPlacement(G4Transform3D(*rotm,
                                                 G4ThreeVector(params->get_double_param("place_x") * cm,
                                                               params->get_double_param("place_y") * cm,
                                                               params->get_double_param("place_z") * cm)),
                                   magnet_logic,
                                   GetName(),
                                   logicMother, false, false, OverlapCheck());

  /* Add volume with solid magnet material */
  G4VSolid *cylinder_solid = new G4Tubs(G4String(GetName().append("_Solid")),
                                        radius,
                                        radius + thickness,
                                        params->get_double_param("length") * cm / 2., 0, twopi);
  G4LogicalVolume *cylinder_logic = new G4LogicalVolume(cylinder_solid,
                                                        TrackerMaterial,
                                                        G4String(GetName()),
                                                        nullptr, nullptr, nullptr);
  cylinder_logic->SetVisAttributes(siliconVis);

  cylinder_physi = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0),
                                     cylinder_logic,
                                     G4String(GetName().append("_Solid")),
                                     magnet_logic, false, false, OverlapCheck());
}
