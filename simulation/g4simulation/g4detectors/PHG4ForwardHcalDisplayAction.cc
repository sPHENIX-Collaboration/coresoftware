#include "PHG4ForwardHcalDisplayAction.h"

#include <g4main/PHG4Utils.h>

#include <Geant4/G4Color.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4VPhysicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>

using namespace std;

PHG4ForwardHcalDisplayAction::PHG4ForwardHcalDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4ForwardHcalDisplayAction::~PHG4ForwardHcalDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4ForwardHcalDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
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
    if (it.second == "Absorber")
    {
      visatt->SetColour(G4Colour::Gray());

    }
    else if (it.second == "FHcalEnvelope")
    {
      visatt->SetVisibility(false);
    }
    else if (it.second == "Scintillator")
    {
      visatt->SetColour(G4Colour::Gray());
    }
    else if (it.second == "SingleScintillator")
    {
      visatt->SetColour(G4Colour::Cyan());
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
