#include "PHG4Detector.h"

#include "PHG4Subsystem.h"

#include <TSystem.h>

#include <Geant4/G4Colour.hh>  // for G4Colour
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4NistManager.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4ThreeVector.hh>     // for G4ThreeVector
#include <Geant4/G4VisAttributes.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/stacktrace.hpp>
#pragma GCC diagnostic pop

PHG4Detector::PHG4Detector(PHG4Subsystem *subsys, PHCompositeNode *Node, const std::string &nam)
  : m_topNode(Node)
  , m_MySubsystem(subsys)
  , m_Name(nam)
{
}

void PHG4Detector::Construct(G4LogicalVolume *world)
{
  PHG4Subsystem *MyMotherSubsystem = m_MySubsystem->GetMotherSubsystem();
  if (MyMotherSubsystem)
  {
    ConstructMe(MyMotherSubsystem->GetLogicalVolume());
  }
  else
  {
    ConstructMe(world);
  }
  return;
}

int PHG4Detector::DisplayVolume(G4VSolid *volume, G4LogicalVolume *logvol, G4RotationMatrix *rotm)
{
  G4LogicalVolume *checksolid = new G4LogicalVolume(volume, G4Material::GetMaterial("G4_POLYSTYRENE"), "DISPLAYLOGICAL", 0, 0, 0);
  int iret = DisplayVolume(checksolid, logvol, rotm);
  return iret;
}

int PHG4Detector::DisplayVolume(G4LogicalVolume *checksolid, G4LogicalVolume *logvol, G4RotationMatrix *rotm)
{
  G4VisAttributes *visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(false);
  switch (m_ColorIndex)
  {
  case 0:
    visattchk->SetColour(G4Colour::Red());
    m_ColorIndex++;
    break;
  case 1:
    visattchk->SetColour(G4Colour::Magenta());
    m_ColorIndex++;
    break;
  case 2:
    visattchk->SetColour(G4Colour::Yellow());
    m_ColorIndex++;
    break;
  case 3:
    visattchk->SetColour(G4Colour::Blue());
    m_ColorIndex++;
    break;
  case 4:
    visattchk->SetColour(G4Colour::Cyan());
    m_ColorIndex++;
    break;
  default:
    visattchk->SetColour(G4Colour::Green());
    m_ColorIndex = 0;
    break;
  }

  checksolid->SetVisAttributes(visattchk);
  new G4PVPlacement(rotm, G4ThreeVector(0, 0, 0), checksolid, "DISPLAYVOL", logvol, 0, false, true);
  return 0;
}

G4Material *PHG4Detector::GetDetectorMaterial(const std::string &name, const bool quit)
{
  G4Material *thismaterial =  G4Material::GetMaterial(name,false);
  if (thismaterial)
  {
    return thismaterial;
  }
  thismaterial = G4NistManager::Instance()->FindOrBuildMaterial(name);
  if (!thismaterial)
  {
    if (!quit)
    {
      return nullptr;
    }
    std::cout << "PHG4Detector::GetDetectorMaterial: Could not locate " << name << " in NIST DB or create it" << std::endl;
    std::cout << boost::stacktrace::stacktrace();
    std::cout << std::endl;
    std::cout << "read the above stack trace who is calling this material" << std::endl;
    gSystem->Exit(1);
    exit(1); // so coverity gets it
  }
  return thismaterial;
}


G4Element *PHG4Detector::GetDetectorElement(const std::string &name, const bool quit)
{
  G4Element *thiselement =  G4Element::GetElement(name,false);
  if (thiselement)
  {
    return thiselement;
  }
  thiselement = G4NistManager::Instance()->FindOrBuildElement(name);
  if (!thiselement)
  {
    if (!quit)
    {
      return nullptr;
    }
    std::cout << "PHG4Detector::GetDetectorElement: Could not locate " << name << " in NIST DB or create it" << std::endl;
    std::cout << boost::stacktrace::stacktrace();
    std::cout << std::endl;
    std::cout << "read the above stack trace who is calling this material" << std::endl;
    gSystem->Exit(1);
    exit(1); // so coverity gets it
  }
  return thiselement;
}
