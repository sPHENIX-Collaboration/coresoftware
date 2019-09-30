#include "G4JLeicDIRCDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>
#include <utility>  // for pair

using namespace std;

G4JLeicDIRCDisplayAction::G4JLeicDIRCDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

G4JLeicDIRCDisplayAction::~G4JLeicDIRCDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void G4JLeicDIRCDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
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
