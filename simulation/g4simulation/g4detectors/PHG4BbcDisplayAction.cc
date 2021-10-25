#include "PHG4BbcDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>
#include <utility>  // for pair

using namespace std;

PHG4BbcDisplayAction::PHG4BbcDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4BbcDisplayAction::~PHG4BbcDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4BbcDisplayAction::ApplyDisplayAction(G4VPhysicalVolume * /*physvol*/)
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
    if (it.second == "Bbc_quartz")
    {
      visatt->SetColour(G4Colour::Yellow());
    }
    else if (it.second == "Bbc_Shell")
    {
      visatt->SetColour(G4Colour::Grey());
    }
    else if (it.second == "Bbc_Front_Plate")
    {
      visatt->SetColour(G4Colour::Grey());
    }
    else
    {
      cout << GetName() << " unknown logical volume " << it.second << endl;
      visatt->SetColour(G4Colour::Cyan());
//      gSystem->Exit(1);
    }
    logvol->SetVisAttributes(visatt);
  }
  return;
}
