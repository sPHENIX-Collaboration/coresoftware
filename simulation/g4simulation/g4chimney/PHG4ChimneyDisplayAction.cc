#include "PHG4ChimneyDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

PHG4ChimneyDisplayAction::PHG4ChimneyDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4ChimneyDisplayAction::~PHG4ChimneyDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4ChimneyDisplayAction::ApplyDisplayAction(G4VPhysicalVolume * /*physvol*/)
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
    m_VisAttVec.push_back(visatt);
    if (it.second == "CHIMNEYBUS" || it.second == "CHIMNEYBUSLEG")
    {
      visatt->SetColour(G4Colour::Brown());
    }
    else if (it.second == "CHIMNEYG10")
    {
      visatt->SetColour(G4Colour::Green());
    }
    else if (it.second == "CHIMNEYDOOR")
    {
      visatt->SetColour(G4Colour::Grey());
    }
    else if (it.second == "CHIMNEYHELINE")
    {
      visatt->SetColour(G4Colour::Blue());
    }
    else if (it.second == "CHIMNEYSHIELD")
    {
      visatt->SetColour(G4Colour::Magenta());
    }
    else if (it.second == "CHIMNEYRISER" || it.second == "CHIMNEYVALVEBOX")
    {
      visatt->SetColour(G4Colour::Yellow());
    }
    else
    {
      visatt->SetColour(G4Colour::Cyan());
    }
    logvol->SetVisAttributes(visatt);
  }
}
