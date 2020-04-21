#include "BeamLineMagnetDetector.h"
#include "BeamLineMagnetDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>  // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>
#include <g4main/PHG4Utils.h>

#include <phool/phool.h>

#include <Geant4/G4ChordFinder.hh>
#include <Geant4/G4ClassicalRK4.hh>
#include <Geant4/G4FieldManager.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Mag_UsualEqRhs.hh>
#include <Geant4/G4MagneticField.hh>
#include <Geant4/G4Material.hh>
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

class G4VSolid;
class PHCompositeNode;
class PHG4Subsystem;

using namespace std;

//_______________________________________________________________
BeamLineMagnetDetector::BeamLineMagnetDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int lyr)
  : PHG4Detector(subsys, Node, dnam)
  , params(parameters)
  , magnet_physi(nullptr)
  , magnet_iron_physi(nullptr)
  , m_DisplayAction(dynamic_cast<BeamLineMagnetDisplayAction *>(subsys->GetDisplayAction()))
  , layer(lyr)
{
}

//_______________________________________________________________
int BeamLineMagnetDetector::IsInBeamLineMagnet(const G4VPhysicalVolume *volume) const
{
  if (volume == magnet_physi)
  {
    return 1;
  }
  else if (volume == magnet_iron_physi)
  {
    return -1;
  }
  return false;
}

//_______________________________________________________________
void BeamLineMagnetDetector::ConstructMe(G4LogicalVolume *logicMother)
{

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
  if (magnettype == "DIPOLE")
  {
    G4ThreeVector field(params->get_double_param("field_x") * tesla,
			params->get_double_param("field_y") * tesla,
			params->get_double_param("field_z") * tesla);
// magnets can be rotated in y
    field.rotateY(params->get_double_param("rot_y") * deg);
    magField = new G4UniformMagField(field);
    if (Verbosity() > 0)
    {
      cout << "Creating DIPOLE with field x: " << field.x() / tesla
	   << ", y: " << field.y() / tesla
	   << ", z: " << field.z()/ tesla
	   << " and name " << GetName() << endl;
    }
  }
  else if (magnettype == "QUADRUPOLE")
  {
    G4double fieldGradient = params->get_double_param("fieldgradient") * tesla / meter;

    /* G4MagneticField::GetFieldValue( pos*, B* ) uses GLOBAL coordinates, not local.
       * Therefore, place magnetic field center at the correct location and angle for the
       * magnet AND do the same transformations for the logical volume (see below). */
    magField = new G4QuadrupoleMagField(fieldGradient, origin, rotm);

    if (Verbosity() > 0)
    {
      cout << "Creating QUADRUPOLE with gradient " << fieldGradient << " and name " << GetName() << endl;
      cout << "at x, y, z = " << origin.x() << " , " << origin.y() << " , " << origin.z() << endl;
      cout << "with rotation around x, y, z axis  by: " << rotm->phiX() << ", " << rotm->phiY() << ", " << rotm->phiZ() << endl;
    }
  }

  if (!magField && Verbosity() > 0)
  {
    cout << PHWHERE << " No magnetic field specified for " << GetName() 
	 << " of type " << magnettype << endl;
  }


  /* Add volume with solid magnet material */
  G4VSolid *magnet_iron_solid = new G4Tubs(G4String(GetName().append("_Solid").c_str()),
                                        0,
					   params->get_double_param("outer_radius") * cm,
                                        params->get_double_param("length") * cm / 2., 0, twopi);
  G4LogicalVolume *magnet_iron_logic = new G4LogicalVolume(magnet_iron_solid,
							G4Material::GetMaterial("G4_Fe"),
                                                        G4String(GetName().c_str()),
                                                        0, 0, 0);
  m_DisplayAction->AddVolume(magnet_iron_logic,magnettype);
  magnet_iron_physi = new G4PVPlacement(G4Transform3D(*rotm,origin),
                                     magnet_iron_logic,
                                     G4String(GetName().append("_Solid").c_str()),
                                     logicMother, 0, false, OverlapCheck());

  G4VSolid *magnet_solid = new G4Tubs(G4String(GetName().c_str()),
                                      0,
                                      params->get_double_param("inner_radius") * cm,
                                      params->get_double_param("length") * cm / 2., 0, twopi);

  G4LogicalVolume *magnet_logic = new G4LogicalVolume(magnet_solid,
                                                      G4Material::GetMaterial("G4_Galactic"),
                                                      G4String(GetName().c_str()),
                                                      0, 0, 0);

  /* Set field manager for logical volume */
  if (magField)
  {
    /* Set up Geant4 field manager */
    G4Mag_UsualEqRhs *localEquation = new G4Mag_UsualEqRhs(magField);
    G4ClassicalRK4 *localStepper = new G4ClassicalRK4(localEquation);
    G4double minStep = 0.25 * mm;  // minimal step, 1 mm is default
    G4ChordFinder *localChordFinder = new G4ChordFinder(magField, minStep, localStepper);

    G4FieldManager *fieldMgr = new G4FieldManager();
    fieldMgr->SetDetectorField(magField);
    fieldMgr->SetChordFinder(localChordFinder);
    magnet_logic->SetFieldManager(fieldMgr, true);
  }
  m_DisplayAction->AddVolume(magnet_logic,"FIELDVOLUME");

  /* create magnet physical volume */

  magnet_physi = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0),
                                   magnet_logic,
                                   G4String(GetName().c_str()),
                                   magnet_iron_logic, 0, false, OverlapCheck());

}
