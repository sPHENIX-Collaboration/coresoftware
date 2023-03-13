#include "PHG4OHCalDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

PHG4OHCalDisplayAction::PHG4OHCalDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4OHCalDisplayAction::~PHG4OHCalDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
  m_ScintiLogVolSet.clear();
}

void PHG4OHCalDisplayAction::ApplyDisplayAction(G4VPhysicalVolume * /*physvol*/)
{
  for (auto &it : m_ScintiLogVolSet)
  {
    if (it->GetVisAttributes())
    {
      return;
    }
    G4VisAttributes *m_VisAtt = new G4VisAttributes();
    m_VisAtt->SetVisibility(true);
    m_VisAtt->SetForceSolid(true);
    m_VisAtt->SetColor(G4Colour::Green());
    it->SetVisAttributes(m_VisAtt);
    m_VisAttVec.push_back(m_VisAtt);
  }

  if (m_SteelVol)
  {
    if (m_SteelVol->GetVisAttributes())
    {
      return;
    }
    G4VisAttributes *m_VisAtt = new G4VisAttributes();
    m_VisAtt->SetVisibility(true);
    m_VisAtt->SetForceSolid(true);
    m_VisAtt->SetColor(G4Colour::Grey());
    m_SteelVol->SetVisAttributes(m_VisAtt);
    m_VisAttVec.push_back(m_VisAtt);
  }

  if (m_ChimSteelVol)
  {
    if (m_ChimSteelVol->GetVisAttributes())
    {
      return;
    }
    G4VisAttributes *m_VisAtt2 = new G4VisAttributes();
    m_VisAtt2->SetVisibility(true);
    m_VisAtt2->SetForceSolid(true);
    m_VisAtt2->SetColor(G4Colour::Grey());
    m_ChimSteelVol->SetVisAttributes(m_VisAtt2);
    m_VisAttVec2.push_back(m_VisAtt2);
  }
}
