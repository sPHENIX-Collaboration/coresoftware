#include "PHG4sPHENIXMagnetDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>
#include <utility>  // for pair

PHG4sPHENIXMagnetDisplayAction::PHG4sPHENIXMagnetDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4sPHENIXMagnetDisplayAction::~PHG4sPHENIXMagnetDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4sPHENIXMagnetDisplayAction::ApplyDisplayAction(G4VPhysicalVolume * /*physvol*/)
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
    if (it.second == "CRYOSTAT" || it.second == "THERM" || it.second == "COILSUP" || it.second == "COIL" || it.second == "CONNECTOR")
    {
      visatt->SetColour(G4Colour::Magenta());
    }
    else if (it.second == "CRYOINT" || it.second == "THERMVAC")
    {
      // visatt->SetVisibility(false);
      // visatt->SetForceSolid(false);
      visatt->SetColour(G4Colour::Gray());
    }
    else if (it.second == "BusbarD")
    {
      visatt->SetColour(G4Colour::Brown());
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
