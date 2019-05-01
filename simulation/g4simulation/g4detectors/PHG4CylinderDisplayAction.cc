#include "PHG4CylinderDisplayAction.h"

#include <g4main/PHG4Utils.h>

#include <phparameter/PHParameters.h>

#include <Geant4/G4VPhysicalVolume.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4VisAttributes.hh>

#include <iostream>

using namespace std;
PHG4CylinderDisplayAction::PHG4CylinderDisplayAction(const std::string &name, PHParameters *pars):
  PHG4DisplayAction(name),
  m_Params(pars),
  m_VisAtt(nullptr)
{}

PHG4CylinderDisplayAction::~PHG4CylinderDisplayAction()
{
  delete m_VisAtt;
}

void PHG4CylinderDisplayAction::ApplyDisplayAction(G4VPhysicalVolume* physvol)
{
// in our case loop over all volumes and see where ours is
  FindVolumes(physvol);

}

int PHG4CylinderDisplayAction::CheckVolume(G4VPhysicalVolume *physvol)
{
  if (physvol == m_MyVolume)
  {
    return CheckReturnCodes::ABORT; //
  }
  return CheckReturnCodes::FAILED;
}

void PHG4CylinderDisplayAction::ApplyVisAttributes(G4VPhysicalVolume *physvol)
{
  G4LogicalVolume* myvol = physvol->GetLogicalVolume();
// check if vis attributes exist, if so someone else has set them and we do nothing
  if (myvol->GetVisAttributes())
  {
    return;
  }
  m_VisAtt = new G4VisAttributes();
  if (m_Params->get_int_param("blackhole"))
  {
    PHG4Utils::SetColour(m_VisAtt, "BlackHole");
    m_VisAtt->SetVisibility(false);
    m_VisAtt->SetForceSolid(false);
  }
  else
  {
    PHG4Utils::SetColour(m_VisAtt, m_Params->get_string_param("material"));
    m_VisAtt->SetVisibility(true);
    m_VisAtt->SetForceSolid(true);
  }
  myvol->SetVisAttributes(m_VisAtt);
  return;
}

