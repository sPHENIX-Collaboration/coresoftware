#include "PHG4TPCDisplayAction.h"

#include <g4main/PHG4ColorDefs.h>

#include <g4main/PHG4Utils.h>

#include <TSystem.h>

#include <Geant4/G4Color.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VPhysicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <iostream>

using namespace std;

PHG4TPCDisplayAction::PHG4TPCDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4TPCDisplayAction::~PHG4TPCDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4TPCDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
{

 static const G4Colour color[] = {PHG4TPCColorDefs::tpc_cu_color,
                                          PHG4TPCColorDefs::tpc_pcb_color,
                                          PHG4TPCColorDefs::tpc_honeycomb_color,
                                          PHG4TPCColorDefs::tpc_cu_color,
                                          PHG4TPCColorDefs::tpc_pcb_color,
                                          PHG4TPCColorDefs::tpc_kapton_color,
                                          PHG4TPCColorDefs::tpc_cu_color,
                                          PHG4TPCColorDefs::tpc_kapton_color,
                                          PHG4TPCColorDefs::tpc_cu_color};
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
    else if (it.second == "TpcEnvelope")
    {
      visatt->SetVisibility(false);
    }
    else if (it.second == "TpcGas")
    {
      visatt->SetColor(PHG4TPCColorDefs::tpc_gas_color);
    }
    else if (it.second == "TpcHoneyComb")
    {
      visatt->SetColor(PHG4TPCColorDefs::tpc_honeycomb_color);
    }
    else if (it.second == "TpcWindow")
    {
      visatt->SetColor(PHG4TPCColorDefs::tpc_pcb_color);
    }
    else
    {
      cout << "did not assing color to " << it.first->GetName()
	   << " under " << it.second << endl;
      gSystem->Exit(1);
    }
    logvol->SetVisAttributes(visatt);
  }
  for (unsigned int i=0; i<m_TpcInnerLayersVec.size(); i++)
  {
    G4VisAttributes *visatt = new G4VisAttributes();
    visatt->SetVisibility(true);
    visatt->SetForceSolid(true);
    m_VisAttVec.push_back(visatt);  // for later deletion
    visatt->SetColor(color[i]);
    m_TpcInnerLayersVec[i]->SetVisAttributes(visatt);
  }
  for (unsigned int i=0; i<m_TpcOuterLayersVec.size(); i++)
  {
    G4VisAttributes *visatt = new G4VisAttributes();
    visatt->SetVisibility(true);
    visatt->SetForceSolid(true);
    m_VisAttVec.push_back(visatt);  // for later deletion
    visatt->SetColor(color[i]);
    m_TpcOuterLayersVec[i]->SetVisAttributes(visatt);
  }

  return;
}
