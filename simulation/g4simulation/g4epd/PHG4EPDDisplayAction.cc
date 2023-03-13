#include "PHG4EPDDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <iostream>  // for operator<<, basic_ostream, endl
#include <string>
#include <utility>  // for pair

PHG4EPDDisplayAction::PHG4EPDDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4EPDDisplayAction::~PHG4EPDDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4EPDDisplayAction::ApplyDisplayAction(G4VPhysicalVolume * /*physvol*/)
{
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
    if (it.second.find("EPD_tile") != std::string::npos)
    {
      visatt->SetColour(G4Colour::Red());
    }
    else
    {
      std::cout << GetName() << " unknown logical volume " << it.second << std::endl;
      visatt->SetColour(G4Colour::White());
    }
    logvol->SetVisAttributes(visatt);
  }
  return;
}
