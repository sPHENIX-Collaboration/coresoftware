#include "PHG4MicromegasDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <iostream>  // for operator<<, basic_ostream, endl
#include <utility>   // for pair
#include <string>

PHG4MicromegasDisplayAction::PHG4MicromegasDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4MicromegasDisplayAction::~PHG4MicromegasDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4MicromegasDisplayAction::ApplyDisplayAction(G4VPhysicalVolume * /*physvol*/)
{
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
    if (it.first->GetName().find("invisible") != std::string::npos)
    {
      visatt->SetColour(it.second);
      visatt->SetVisibility(false);
    }
    else
    {
      visatt->SetColour(it.second);
    }
    logvol->SetVisAttributes(visatt);
  }
  return;
}
