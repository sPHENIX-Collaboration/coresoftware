#include "PHG4TpcEndCapDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Utils.h>

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4VisAttributes.hh>

#include <utility>  // for pair

PHG4TpcEndCapDisplayAction::PHG4TpcEndCapDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4TpcEndCapDisplayAction::~PHG4TpcEndCapDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4TpcEndCapDisplayAction::ApplyDisplayAction(G4VPhysicalVolume * /*physvol*/)
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

    if (it.second == "G10" or it.second == "FR4")
    {
      visatt->SetColour(G4Colour(0.0, .5, 0.0));
    }
    else if (it.second == "wagon_wheel")
    {
      visatt->SetColour(G4Colour(.8, 0.0, 0.0));
      visatt->SetForceLineSegmentsPerCircle(100);
    }
    else if (it.second == "cooling_block")
    {
      visatt->SetColour(G4Colour(.8, .8, .8));
      visatt->SetForceLineSegmentsPerCircle(100);
    }
    else
    {
      PHG4Utils::SetColour(visatt, it.first->GetMaterial()->GetName());
    }
    logvol->SetVisAttributes(visatt);
  }
  return;
}
