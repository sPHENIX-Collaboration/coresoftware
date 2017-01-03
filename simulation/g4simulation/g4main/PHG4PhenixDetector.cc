#include "PHG4PhenixDetector.h"
#include "PHG4Detector.h"
#include "PHG4RegionInformation.h"

#include <phool/recoConsts.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Element.hh>
#include <Geant4/G4GeometryManager.hh>
#include <Geant4/G4LogicalVolumeStore.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PhysicalVolumeStore.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4Region.hh>
#include <Geant4/G4RegionStore.hh>
#include <Geant4/G4SolidStore.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <iostream>

using namespace std;

//____________________________________________________________________________
PHG4PhenixDetector::PHG4PhenixDetector( void ):
  verbosity(0),
  defaultMaterial( NULL ),
  logicWorld( NULL ),
  physiWorld( NULL ),
  WorldSizeX(1000*cm),
  WorldSizeY(1000*cm),
  WorldSizeZ(1000*cm),
  worldshape("G4TUBS"),
  worldmaterial("G4_AIR")
{
}

PHG4PhenixDetector::~PHG4PhenixDetector()
{
  while (detectors_.begin() != detectors_.end())
    {
      delete detectors_.back();
      detectors_.pop_back();
    }
}


//_______________________________________________________________________________________________
G4VPhysicalVolume* PHG4PhenixDetector::Construct()
{

  recoConsts *rc = recoConsts::instance();
  if (verbosity > 0) std::cout << "PHG4PhenixDetector::Construct." << std::endl;
  // Clean old geometry, if any
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  if (verbosity > 0) std::cout << "PHG4PhenixDetector::Construct - cleaning done." << std::endl;
  //default materials of the World
  //  defaultMaterial  = nist->FindOrBuildMaterial("G4_AIR");

  // World
  G4VSolid *solidWorld = NULL;
  if (worldshape == "G4BOX")
    {
      solidWorld = new G4Box("World", WorldSizeX / 2, WorldSizeY / 2, WorldSizeZ / 2);
    }
  else if (worldshape == "G4Tubs")
    {
      solidWorld = new G4Tubs("World", 0., WorldSizeY / 2, WorldSizeZ / 2, 0, 2 * M_PI);
    }
  else
    {
      cout << "Unknown world shape " << worldshape << endl;
      cout << "implemented are G4BOX, G4Tubs" << endl;
      exit(1);
    }
  rc->set_CharFlag("WorldShape",solidWorld->GetEntityType()); // needed for checks if a particle is inside or outside of our world
  logicWorld = new G4LogicalVolume(solidWorld, G4Material::GetMaterial(worldmaterial),"World");
  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
  physiWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0 );

  G4Region* defaultRegion = (*(G4RegionStore::GetInstance()))[0];
  PHG4RegionInformation* info = new PHG4RegionInformation();
  info->SetWorld();
  defaultRegion->SetUserInformation(info);
  if (verbosity > 0) {
    std::cout << "PHG4PhenixDetector::Construct " << solidWorld->GetEntityType() << " world "
	      << "material " << logicWorld->GetMaterial()->GetName() << " done." << std::endl;
  }
  
  // construct all detectors
  for( DetectorList::iterator iter = detectors_.begin(); iter != detectors_.end(); ++iter )
  {
    if( *iter )
    {
      (*iter)->Construct( logicWorld );
    }
  }

  if (verbosity > 0) std::cout << "PHG4PhenixDetector::Construct - done." << std::endl;

  return physiWorld;
}

