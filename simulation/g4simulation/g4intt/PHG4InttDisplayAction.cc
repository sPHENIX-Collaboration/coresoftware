#include "PHG4InttDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <TSystem.h>

#include <Geant4/G4Colour.hh>          // for G4Colour
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4VisAttributes.hh>

#include <iostream>
#include <utility>                     // for pair

using namespace std;

PHG4InttDisplayAction::PHG4InttDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4InttDisplayAction::~PHG4InttDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4InttDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
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
    if (it.second == "FPHX")
    {
      visatt->SetColour(G4Colour::Blue());
    }
    else if (it.second == "HdiCopper")
    {
      visatt->SetColour(G4Colour::White());
    }
    else if (it.second == "HdiKapton")
    {
      visatt->SetColour(G4Colour::Yellow());
    }
    else if (it.second == "Ladder")
    {
      visatt->SetVisibility(false);
      visatt->SetForceSolid(false);
    }
    else if (it.second == "PGS")
    {
      visatt->SetColour(G4Colour::Red());
    }
    else if (it.second == "Rail")
    {
      visatt->SetColour(G4Colour::Cyan());
    }
    else if (it.second == "RohaCell")
    {
      visatt->SetColour(G4Colour::White());
    }
    else if (it.second == "SiActive")
    {
      visatt->SetColour(G4Colour::White());
    }
    else if (it.second == "SiInActive")
    {
      visatt->SetColour(G4Colour::Red());
    }
    else if (it.second == "StaveBox")
    {
      visatt->SetVisibility(false);
      visatt->SetForceSolid(false);
    }
    else if (it.second == "StaveCooler")
    {
      visatt->SetColour(G4Colour::Grey());
    }
    else if (it.second == "StaveCurve")
    {
      visatt->SetColour(G4Colour::Grey());
    }
    else if (it.second == "StaveGlueBox")
    {
      visatt->SetColour(G4Colour::Cyan());
    }
    else if (it.second == "StavePipe")
    {
      visatt->SetColour(G4Colour::White());
    }
    else if (it.second == "StaveStraightInner")
    {
      visatt->SetColour(G4Colour::Grey());
    }
    else if (it.second == "StaveStraightOuter")
    {
      visatt->SetColour(G4Colour::Grey());
    }
    else if (it.second == "StaveWater")
    {
      visatt->SetColour(G4Colour::Blue());
    }
    else
    {
      cout << "did not assing color to " << it.first->GetName()
	   << " under " << it.second << endl;
      gSystem->Exit(1);
    }
    logvol->SetVisAttributes(visatt);
  }
  return;
}
