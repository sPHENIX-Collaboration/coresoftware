#include "PHG4BbcDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#include <iostream>
#include <utility>  // for pair

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
    if (it.second == "Bbc_Breeder_Module")
    {
      visatt->SetColour(G4Colour::Green());
    }
    else if (it.second == "Bbc_Cover_Plates")
    {
      visatt->SetColour(G4Colour::Gray());
    }
    else if (it.second == "Bbc_Inner_Shell")
    {
      visatt->SetColour(G4Colour::Gray());
    }
    else if (it.second == "Bbc_quartz")
    {
      visatt->SetColour(G4Colour::Cyan());
    }
    else if (it.second == "Bbc_Outer_Shell")
    {
      visatt->SetColour(G4Colour::Gray());
    }
    else if (it.second == "Bbc_PMT")
    {
      visatt->SetColour(G4Colour::Blue());
    }
    else if (it.second == "Bbc_Shell")
    {
      visatt->SetColour(0.5, 0.5, 0.5, 0.4);
    }
    else if (it.second == "Bbc_Front_Plate")
    {
      visatt->SetColour(0.5, 0.5, 0.5, 0.4);
    }
    else if (it.second == "Bbc_attach_plate")
    {
      visatt->SetColour(G4Colour::Gray());
    }
    else if (it.second == "Bbc_CableCond")
    {
      visatt->SetColour(G4Colour::Yellow());
    }
    else if (it.second == "Bbc_CableShield")
    {
      visatt->SetColour(G4Colour::White());
    }
    else if (it.second == "Bbc_Base_Plates")
    {
      visatt->SetColour(G4Colour::Gray());
    }
    else
    {
      std::cout << "PHG4BbcDisplayAction::ApplyDisplayAction unknown logical volume " << it.second << " in " << GetName() << std::endl;
      gSystem->Exit(1);
    }
    logvol->SetVisAttributes(visatt);
  }
  return;
}
