#include "PHG4RICHDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>
#include <utility>                     // for pair

using namespace std;

PHG4RICHDisplayAction::PHG4RICHDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4RICHDisplayAction::~PHG4RICHDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4RICHDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
{
  // check if vis attributes exist, if so someone else has set them and we do nothing
  for (auto it : m_LogicalVolumeMap)
  {
    G4LogicalVolume *logvol = it.first;
    if (logvol->GetVisAttributes())
    {
      continue;
    }
    G4VisAttributes *visatt = new G4VisAttributes();
    visatt->SetVisibility(true);
    visatt->SetForceSolid(true);
    m_VisAttVec.push_back(visatt);  // for later deletion
    if (it.second == "HBD")
    {
      visatt->SetColour(G4Colour::Red());
    }
    else if (it.second == "Mirror")
    {
      visatt->SetColour(G4Colour::Green());
      visatt->SetForceLineSegmentsPerCircle(50);
    }
    else if (it.second == "Sector")
    {
      visatt->SetColour(G4Colour::White());
      visatt->SetForceWireframe(true);
      visatt->SetForceLineSegmentsPerCircle(50);
    }
    else if (it.second == "Window")
    {
      visatt->SetColour(G4Colour::Yellow());
      visatt->SetForceLineSegmentsPerCircle(50);
    }
    else
    {
      cout << "unknown logical volume " << it.second << endl;
      gSystem->Exit(1);
    }
    logvol->SetVisAttributes(visatt);
  }
  return;
}
