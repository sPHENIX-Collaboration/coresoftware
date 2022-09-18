#include "BeamLineMagnetDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <phool/phool.h>

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>
#include <utility>  // for pair

BeamLineMagnetDisplayAction::BeamLineMagnetDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

BeamLineMagnetDisplayAction::~BeamLineMagnetDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void BeamLineMagnetDisplayAction::ApplyDisplayAction(G4VPhysicalVolume * /*physvol*/)
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
    if (it.second == "DIPOLE")
    {
      visatt->SetColour(G4Color(0.2, 0.8, 0.2, 0.8));
    }
    else if (it.second == "FIELDVOLUME" || it.second == "OFF")
    {
      visatt->SetVisibility(false);
      visatt->SetForceSolid(false);
    }
    else if (it.second == "QUADRUPOLE")
    {
      visatt->SetColour(G4Color(0., 0.3, 0.7, 0.8));
    }
    else if (it.second == "SEXTUPOLE")
    {
      visatt->SetColour(G4Color(1., 0.5, 0.7, 0.8));
    }
    else
    {
      std::cout << PHWHERE << "unknown logical volume " << it.second << std::endl;
      gSystem->Exit(1);
    }
    logvol->SetVisAttributes(visatt);
  }
  return;
}
