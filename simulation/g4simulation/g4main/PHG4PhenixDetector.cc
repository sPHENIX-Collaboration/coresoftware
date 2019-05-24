#include "PHG4PhenixDetector.h"
#include "PHG4Detector.h"
#include "PHG4PhenixDisplayAction.h"
#include "PHG4Reco.h"
#include "PHG4RegionInformation.h"

#include <phool/recoConsts.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Element.hh>
#include <Geant4/G4GeometryManager.hh>
#include <Geant4/G4LogicalVolumeStore.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalVolumeStore.hh>
#include <Geant4/G4Region.hh>
#include <Geant4/G4RegionStore.hh>
#include <Geant4/G4SolidStore.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Tubs.hh>

#include <boost/foreach.hpp>

#include <cmath>
#include <iostream>

using namespace std;

//____________________________________________________________________________
PHG4PhenixDetector::PHG4PhenixDetector(PHG4Reco *subsys)
  : m_DisplayAction(dynamic_cast<PHG4PhenixDisplayAction *>(subsys->GetDisplayAction()))
  , m_Verbosity(0)
  , defaultMaterial(nullptr)
  , logicWorld(nullptr)
  , physiWorld(nullptr)
  , WorldSizeX(1000 * cm)
  , WorldSizeY(1000 * cm)
  , WorldSizeZ(1000 * cm)
  , worldshape("G4TUBS")
  , worldmaterial("G4_AIR")
{
}

PHG4PhenixDetector::~PHG4PhenixDetector()
{
  while (m_DetectorList.begin() != m_DetectorList.end())
  {
    delete m_DetectorList.back();
    m_DetectorList.pop_back();
  }
}

//_______________________________________________________________________________________________
G4VPhysicalVolume *PHG4PhenixDetector::Construct()
{
  recoConsts *rc = recoConsts::instance();
  if (m_Verbosity > 0) std::cout << "PHG4PhenixDetector::Construct." << std::endl;
  // Clean old geometry, if any
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  if (m_Verbosity > 0) std::cout << "PHG4PhenixDetector::Construct - cleaning done." << std::endl;
  //default materials of the World
  //  defaultMaterial  = nist->FindOrBuildMaterial("G4_AIR");

  // World
  G4VSolid *solidWorld = nullptr;
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
  rc->set_CharFlag("WorldShape", solidWorld->GetEntityType());  // needed for checks if a particle is inside or outside of our world
  logicWorld = new G4LogicalVolume(solidWorld, G4Material::GetMaterial(worldmaterial), "World");
  m_DisplayAction->AddVolume(logicWorld, "World");
  physiWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0);

  G4Region *defaultRegion = (*(G4RegionStore::GetInstance()))[0];
  PHG4RegionInformation *info = new PHG4RegionInformation();
  info->SetWorld();
  defaultRegion->SetUserInformation(info);
  if (m_Verbosity > 0)
  {
    std::cout << "PHG4PhenixDetector::Construct " << solidWorld->GetEntityType() << " world "
              << "material " << logicWorld->GetMaterial()->GetName() << " done." << std::endl;
  }

  // construct all detectors
  BOOST_FOREACH (PHG4Detector *det, m_DetectorList)
  {
    if (det)
    {
      det->Construct(logicWorld);
    }
  }

  if (m_Verbosity > 0) std::cout << "PHG4PhenixDetector::Construct - done." << std::endl;

  return physiWorld;
}
