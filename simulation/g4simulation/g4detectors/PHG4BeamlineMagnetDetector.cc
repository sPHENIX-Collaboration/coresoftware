#include "PHG4BeamlineMagnetDetector.h"
#include "PHG4Parameters.h"

#include <g4main/PHG4PhenixDetector.h>
#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/phool.h>

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4FieldManager.hh>
#include <Geant4/G4MagneticField.hh>
#include <Geant4/G4UniformMagField.hh>
#include <Geant4/G4QuadrupoleMagField.hh>
#include <Geant4/G4Mag_UsualEqRhs.hh>
#include <Geant4/G4ChordFinder.hh>
#include <Geant4/G4ClassicalRK4.hh>
#include <Geant4/G4TransportationManager.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4VisAttributes.hh>


using namespace std;

//_______________________________________________________________
PHG4BeamlineMagnetDetector::PHG4BeamlineMagnetDetector( PHCompositeNode *Node,  PHG4Parameters *parameters, const std::string &dnam, const int lyr ):
  PHG4Detector(Node,dnam),
  params(parameters),
  magnet_physi(NULL),
  cylinder_physi(NULL),
  layer(lyr)
{}

//_______________________________________________________________
bool PHG4BeamlineMagnetDetector::IsInBeamlineMagnet(const G4VPhysicalVolume * volume) const
{
  if (volume == magnet_physi)
    {
      return true;
    }
  return false;
}


//_______________________________________________________________
void PHG4BeamlineMagnetDetector::Construct( G4LogicalVolume* logicWorld )
{
  G4Material *TrackerMaterial = G4Material::GetMaterial(params->get_string_param("material"));

  if ( ! TrackerMaterial)
    {
      std::cout << "Error: Can not set material" << std::endl;
      exit(-1);
    }

  G4VisAttributes* fieldVis= new G4VisAttributes();
  PHG4Utils::SetColour(fieldVis, "BlackHole");
  fieldVis->SetVisibility(false);
  fieldVis->SetForceSolid(false);

  G4VisAttributes* siliconVis= new G4VisAttributes();
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
  G4ThreeVector origin(params->get_double_param("place_x")*cm,
                       params->get_double_param("place_y")*cm,
                       params->get_double_param("place_z")*cm);

  /* Define magnet rotation matrix */
  G4RotationMatrix *rotm = new G4RotationMatrix();
  rotm->rotateX(params->get_double_param("rot_x")*deg);
  rotm->rotateY(params->get_double_param("rot_y")*deg);
  rotm->rotateZ(params->get_double_param("rot_z")*deg);

  /* Creating a magnetic field */
  G4MagneticField* magField = NULL;

  string magnettype = params->get_string_param("magtype");
  if ( magnettype == "dipole" )
    {
      G4double fieldValue = params->get_double_param("field_y")*tesla;
      magField = new G4UniformMagField(G4ThreeVector(0.,fieldValue,0.));

      if ( verbosity > 0 )
        cout << "Creating DIPOLE with field " << fieldValue << " and name " << GetName() << endl;
    }
  else if ( magnettype == "quadrupole" )
    {
      G4double fieldGradient = params->get_double_param("fieldgradient")*tesla/meter;

      /* G4MagneticField::GetFieldValue( pos*, B* ) uses GLOBAL coordinates, not local.
       * Therefore, place magnetic field center at the correct location and angle for the
       * magnet AND do the same transformations for the logical volume (see below). */
      magField = new G4QuadrupoleMagField ( fieldGradient, origin, rotm );
      //      magField = new PHG4QuadrupoleMagField ( fieldGradient, origin, rotm );

      if ( verbosity > 0 )
        {
          cout << "Creating QUADRUPOLE with gradient " << fieldGradient << " and name " << GetName() << endl;
          cout << "at x, y, z = " << origin.x() << " , " << origin.y() << " , " << origin.z() << endl;
          cout << "with rotation around x, y, z axis  by: " << rotm->phiX() << ", " << rotm->phiY() << ", " << rotm->phiZ() << endl;
        }
    }

  if ( !magField )
    {
      cout << PHWHERE << " Aborting, no magnetic field specified for " << GetName() << endl;
      exit(1);
    }

  /* Set up Geant4 field manager */
  G4Mag_UsualEqRhs* localEquation = new G4Mag_UsualEqRhs(magField);
  G4ClassicalRK4* localStepper = new G4ClassicalRK4( localEquation );
  G4double minStep = 0.25*mm; // minimal step, 1 mm is default
  G4ChordFinder* localChordFinder = new G4ChordFinder( magField, minStep, localStepper );

  G4FieldManager* fieldMgr = new G4FieldManager();
  fieldMgr->SetDetectorField(magField);
  fieldMgr->SetChordFinder( localChordFinder );

  /* Add volume with magnetic field */
  double radius = params->get_double_param("radius")*cm;
  double thickness = params->get_double_param("thickness")*cm;

  G4VSolid *magnet_solid = new G4Tubs(G4String(GetName().c_str()),
                                      0,
                                      radius+thickness,
                                      params->get_double_param("length")*cm/2. ,0,twopi);

  G4LogicalVolume *magnet_logic = new G4LogicalVolume(magnet_solid,
                                                      G4Material::GetMaterial("G4_Galactic"),
                                                      G4String(GetName().c_str()),
                                                      0,0,0);
  magnet_logic->SetVisAttributes(fieldVis);

  /* Set field manager for logical volume */
  G4bool allLocal = true;
  magnet_logic->SetFieldManager(fieldMgr,allLocal);



  /* create magnet physical volume */
  magnet_physi = new G4PVPlacement(G4Transform3D(*rotm,
                                                 G4ThreeVector(params->get_double_param("place_x")*cm,
                                                               params->get_double_param("place_y")*cm,
                                                               params->get_double_param("place_z")*cm) ),
                                   magnet_logic,
                                   G4String(GetName().c_str()),
                                   logicWorld, 0, false, overlapcheck);


  /* Add volume with solid magnet material */
  G4VSolid *cylinder_solid = new G4Tubs(G4String(GetName().append("_Solid").c_str()),
                                        radius,
                                        radius + thickness,
                                        params->get_double_param("length")*cm/2. ,0,twopi);
  G4LogicalVolume *cylinder_logic = new G4LogicalVolume(cylinder_solid,
                                                        TrackerMaterial,
                                                        G4String(GetName().c_str()),
                                                        0,0,0);
  cylinder_logic->SetVisAttributes(siliconVis);

  cylinder_physi = new G4PVPlacement(0, G4ThreeVector(0,0,0),
                                     cylinder_logic,
                                     G4String(GetName().append("_Solid").c_str()),
                                     magnet_logic, 0, false, overlapcheck);

}
