#include "PHG4PhenixDisplayAction.h"

#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>
#include <utility>                     // for pair

using namespace std;

PHG4PhenixDisplayAction::PHG4PhenixDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4PhenixDisplayAction::~PHG4PhenixDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4PhenixDisplayAction::ApplyDisplayAction(G4VPhysicalVolume */*physvol*/)
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
    m_VisAttVec.push_back(visatt);  // for later deletion
    if (it.second == "World")
    {
      visatt->SetVisibility(false);
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
