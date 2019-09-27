#include "PHG4Detector.h"

#include "PHG4Subsystem.h"

#include <Geant4/G4Colour.hh>            // for G4Colour
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4ThreeVector.hh>       // for G4ThreeVector
#include <Geant4/G4VisAttributes.hh>

PHG4Detector::PHG4Detector(PHG4Subsystem *subsys, PHCompositeNode *Node, const std::string &nam)
  : m_topNode(Node)
  , m_MySubsystem(subsys)
  , m_Verbosity(0)
  , m_OverlapCheck(false)
  , m_ColorIndex(0)
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
