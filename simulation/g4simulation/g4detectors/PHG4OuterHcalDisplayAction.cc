#include "PHG4OuterHcalDisplayAction.h"

#include <Geant4/G4Color.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VPhysicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <iostream>

using namespace std;

PHG4OuterHcalDisplayAction::PHG4OuterHcalDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
  , m_MyTopVolume(nullptr)
  , m_SteelVol(nullptr)
{
}

PHG4OuterHcalDisplayAction::~PHG4OuterHcalDisplayAction()
{
  for (auto &it:m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
  m_ScintiLogVolSet.clear();
}

void PHG4OuterHcalDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
{
  for (auto &it:m_ScintiLogVolSet)
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
