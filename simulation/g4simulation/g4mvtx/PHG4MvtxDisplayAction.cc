#include "PHG4MvtxDisplayAction.h"

#include "g4main/PHG4DisplayAction.h"  // for PHG4DisplayAction
#include <g4main/PHG4Utils.h>

#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <utility>                     // for pair

using namespace std;

PHG4MvtxDisplayAction::PHG4MvtxDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4MvtxDisplayAction::~PHG4MvtxDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4MvtxDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *)
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
    if (it.second == "Carbon")
    {
      visatt->SetColour(0.5, 0.5, 0.5, .25);
    }
    else if (it.second == "M60J3K")
    {
      visatt->SetColour(0.25, 0.25, 0.25, .25);
    }
    else if (it.second == "WATER")
    {
      visatt->SetColour(0.0, 0.5, 0.0, .25);
    }
    else if (it.second == "SI")
    {
      PHG4Utils::SetColour(visatt, "G4_Si");
    }
    else if (it.second == "KAPTON")
    {
      PHG4Utils::SetColour(visatt, "G4_KAPTON");
    }
    else if (it.second == "ALUMINUM")
    {
      PHG4Utils::SetColour(visatt, "G4_Al");
    }
    else
    {
      visatt->SetColour(.2, .2, .7, .25);
    }
    visatt->SetVisibility(true);
    visatt->SetForceSolid(true);
    logvol->SetVisAttributes(visatt);
  }
  return;
}
