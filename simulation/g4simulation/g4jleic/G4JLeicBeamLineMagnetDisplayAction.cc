#include "G4JLeicBeamLineMagnetDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <phool/phool.h>

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>
#include <utility>  // for pair

using namespace std;

G4JLeicBeamLineMagnetDisplayAction::G4JLeicBeamLineMagnetDisplayAction(const std::string &name, PHParameters *pars)
  : PHG4DisplayAction(name)
  , m_Params(pars)
{
}

G4JLeicBeamLineMagnetDisplayAction::~G4JLeicBeamLineMagnetDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void G4JLeicBeamLineMagnetDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
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
    if (it.second == "QUADRUPOLE")
    {
      visatt->SetColour(G4Color(0.8, 0.3, 0.1, 0.9));
    }
    else if (it.second == "FIELDVOLUME")
    {
      visatt->SetColour(G4Colour::Gray());
      visatt->SetForceSolid(false);
    }
    else if (it.second == "DIPOLE")
    {
      visatt->SetColour(G4Color(0.2, 0.8, 0.2, 1.));
    }
    else if (it.second == "SOLENOID")
    {
      visatt->SetColour(G4Color(1., 0.5, 0.7, 1.));
    }
    else
    {
      cout << PHWHERE << "unknown logical volume " << it.second << endl;
      gSystem->Exit(1);
    }
    logvol->SetVisAttributes(visatt);
  }
  return;
}
