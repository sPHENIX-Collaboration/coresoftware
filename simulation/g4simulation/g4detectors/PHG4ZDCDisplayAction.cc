#include "PHG4ZDCDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>
#include <utility>  // for pair

PHG4ZDCDisplayAction::PHG4ZDCDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4ZDCDisplayAction::~PHG4ZDCDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4ZDCDisplayAction::ApplyDisplayAction(G4VPhysicalVolume * /*physvol*/)
{
  // check if vis attributes exist, if so someone else has set them and we do nothing
  for (const auto &it : m_LogicalVolumeMap)
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
    if (it.second == "Absorber")
    {
      visatt->SetColour(G4Colour::Blue());
    }
    else if (it.second == "Envelope" || it.second == "fiber_plate_air")
    {
      visatt->SetVisibility(false);
      visatt->SetForceSolid(false);
    }
    else if (it.second == "Fiber")
    {
      visatt->SetColour(G4Colour::Cyan());
    }
    else if (it.second == "FrontBackPlate")
    {
      visatt->SetColour(G4Colour::Red());
    }
    else if (it.second == "Window")
    {
      visatt->SetColour(G4Colour::Blue());
    }
    else if (it.second == "SMD")
    {
      visatt->SetColour(G4Colour::Yellow());
    }
    else if (it.second == "FiberPlate")
    {
      visatt->SetColour(G4Colour::Cyan());
    }
    else if (it.second == "Scint_solid")
    {
      visatt->SetColour(G4Colour::Cyan());
    }
    else
    {
      std::cout << GetName() << " unknown logical volume " << it.second << std::endl;
      gSystem->Exit(1);
    }
    logvol->SetVisAttributes(visatt);
  }
  return;
}
