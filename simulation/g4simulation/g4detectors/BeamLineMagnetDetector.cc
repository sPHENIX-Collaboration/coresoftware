#include "BeamLineMagnetDetector.h"
#include "BeamLineMagnetDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

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

#include <CLHEP/Units/SystemOfUnits.h>  // for cm, deg, tesla, twopi, meter

#include <cassert>
#include <iostream>  // for operator<<, basic_ostream

class G4VSolid;
class PHCompositeNode;

using namespace std;

//_______________________________________________________________
BeamLineMagnetDetector::BeamLineMagnetDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int magnet_id)
  : PHG4Detector(subsys, Node, dnam)
  , params(parameters)
  , magnet_physi(nullptr)
  , magnet_iron_physi(nullptr)
  , m_DisplayAction(dynamic_cast<BeamLineMagnetDisplayAction *>(subsys->GetDisplayAction()))
  , m_MagnetId(magnet_id)
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

  // creating mother volume (cylinder for the time being)
  G4VSolid *magnet_mother_solid = new G4Tubs(G4String(GetName().append("_Mother").c_str()),
                                             0,
                                             params->get_double_param("outer_radius") * cm,
                                             params->get_double_param("length") * cm / 2., 0, twopi);
  G4LogicalVolume *magnet_mother_logic = new G4LogicalVolume(magnet_mother_solid,
                                                             G4Material::GetMaterial("G4_Galactic"),
                                                             G4String(GetName().c_str()),
                                                             0, 0, 0);

  m_DisplayAction->AddVolume(magnet_mother_logic, "OFF");
  new G4PVPlacement(G4Transform3D(*rotm, origin),
                    magnet_mother_logic,
                    G4String(GetName().append("_Mother").c_str()),
                    logicMother, false, m_MagnetId, OverlapCheck());

  G4ThreeVector field_origin(origin);
  G4RotationMatrix *field_rotm = rotm;

  // mother subsys ! = the world and therefore need to explicitly assign global geometry for field mamnagers
  if (GetMySubsystem()->GetMotherSubsystem())
  {
    if (Verbosity() > 0)
    {
      cout << __PRETTY_FUNCTION__ << ": set field using the global coordinate system, as the magnet is a daughter vol. of "
           << GetMySubsystem()->GetMotherSubsystem()->Name() << endl;
    }
    /* Define origin vector (center of magnet) */
    // abs. position to world for field manager
    field_origin = G4ThreeVector(params->get_double_param("field_global_position_x") * cm,
                                 params->get_double_param("field_global_position_y") * cm,
                                 params->get_double_param("field_global_position_z") * cm);

    /* Define magnet rotation matrix */
    // abs. rotation to world for field manager
    field_rotm = new G4RotationMatrix();
    rotm->rotateX(params->get_double_param("field_global_rot_x") * deg);
    rotm->rotateY(params->get_double_param("field_global_rot_y") * deg);
    rotm->rotateZ(params->get_double_param("field_global_rot_z") * deg);
  }

  string magnettype = params->get_string_param("magtype");
  if (magnettype == "DIPOLE")
  {
    G4ThreeVector field(params->get_double_param("field_x") * tesla,
                        params->get_double_param("field_y") * tesla,
                        params->get_double_param("field_z") * tesla);
    // magnets can be rotated in y
    //     field.rotateY(params->get_double_param("rot_y") * deg);
    // apply consistent geometry ops to field and detectors
    field.transform(*field_rotm);
    m_magField = new G4UniformMagField(field);
    if (Verbosity() > 0)
    {
      cout << "Creating DIPOLE with field x: " << field.x() / tesla
           << ", y: " << field.y() / tesla
           << ", z: " << field.z() / tesla
           << " and name " << GetName() << endl;
    }
  }
  else if (magnettype == "QUADRUPOLE")
  {
    G4double fieldGradient = params->get_double_param("fieldgradient") * tesla / meter;

    /* G4MagneticField::GetFieldValue( pos*, B* ) uses GLOBAL coordinates, not local.
       * Therefore, place magnetic field center at the correct location and angle for the
       * magnet AND do the same transformations for the logical volume (see below). */
    m_magField = new G4QuadrupoleMagField(fieldGradient, field_origin, field_rotm);

    if (Verbosity() > 0)
    {
      cout << "Creating QUADRUPOLE with gradient " << fieldGradient << " and name " << GetName() << endl;
      cout << "at x, y, z = " << field_origin.x() << " , " << field_origin.y() << " , " << field_origin.z() << endl;
      cout << "with rotation around x, y, z axis  by: " << field_rotm->phiX() << ", " << field_rotm->phiY() << ", " << field_rotm->phiZ() << endl;
    }
  }
  else
  {
    cout << __PRETTY_FUNCTION__ << " : Fatal error, magnet type of " << magnettype << " is not yet supported" << endl;
    exit(1);
  }

  if (!m_magField && Verbosity() > 0)
  {
    cout << PHWHERE << " No magnetic field specified for " << GetName()
         << " of type " << magnettype << endl;
  }

  /* Add volume with solid magnet material */
  G4VSolid *magnet_iron_solid = new G4Tubs(G4String(GetName().append("_Solid").c_str()),
                                           params->get_double_param("inner_radius") * cm,
                                           params->get_double_param("outer_radius") * cm,
                                           params->get_double_param("length") * cm / 2., 0, twopi);
  G4LogicalVolume *magnet_iron_logic = new G4LogicalVolume(magnet_iron_solid,
                                                           G4Material::GetMaterial("G4_Fe"),
                                                           G4String(GetName().c_str()),
                                                           0, 0, 0);
  m_DisplayAction->AddVolume(magnet_iron_logic, magnettype);
  magnet_iron_physi = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0),
                                        magnet_iron_logic,
                                        G4String(GetName().append("_Solid").c_str()),
                                        magnet_mother_logic, false, m_MagnetId, OverlapCheck());

  G4VSolid *magnet_field_solid = new G4Tubs(G4String(GetName().append("_Field_Solid").c_str()),
                                            0,
                                            params->get_double_param("inner_radius") * cm,
                                            params->get_double_param("length") * cm / 2., 0, twopi);

  m_magnetFieldLogic = new G4LogicalVolume(magnet_field_solid,
                                           G4Material::GetMaterial("G4_Galactic"),
                                           G4String(GetName().c_str()),
                                           0, 0, 0);

  // allow installing new detector inside the magnet, magnet_field_logic
  PHG4Subsystem *mysys = GetMySubsystem();
  mysys->SetLogicalVolume(m_magnetFieldLogic);

  m_DisplayAction->AddVolume(m_magnetFieldLogic, "FIELDVOLUME");

  /* create magnet physical volume */

  magnet_physi = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0),
                                   m_magnetFieldLogic,
                                   G4String(GetName().c_str()),
                                   magnet_mother_logic, false, m_MagnetId, OverlapCheck());
}

//! Optional PostConstruction call after all geometry is constructed
void BeamLineMagnetDetector::PostConstruction()
{
  /* Set field manager for logical volume */
  // This has to be in PostConstruction() in order to apply to all daughter vol. including daughter subsystems

  assert(m_magField);
  assert(m_magnetFieldLogic);

  if (Verbosity() > 0)
  {
    cout << __PRETTY_FUNCTION__ << ": set field to vol " << m_magnetFieldLogic->GetName() << " that include "
         << m_magnetFieldLogic->GetNoDaughters() << " daughter vols." << endl;
  }

  /* Set up Geant4 field manager */
  G4Mag_UsualEqRhs *localEquation = new G4Mag_UsualEqRhs(m_magField);
  G4ClassicalRK4 *localStepper = new G4ClassicalRK4(localEquation);
  G4double minStep = 0.25 * mm;  // minimal step, 1 mm is default
  G4ChordFinder *localChordFinder = new G4ChordFinder(m_magField, minStep, localStepper);

  G4FieldManager *fieldMgr = new G4FieldManager();
  fieldMgr->SetDetectorField(m_magField);
  fieldMgr->SetChordFinder(localChordFinder);
  m_magnetFieldLogic->SetFieldManager(fieldMgr, true);
}
