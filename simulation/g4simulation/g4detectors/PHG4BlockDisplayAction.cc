#include "PHG4BlockDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>
#include <g4main/PHG4Utils.h>

#include <phparameter/PHParameters.h>

#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

using namespace std;

PHG4BlockDisplayAction::PHG4BlockDisplayAction(const std::string &name, PHParameters *pars)
  : PHG4DisplayAction(name)
  , m_Params(pars)
  , m_VisAtt(nullptr)
{
}

PHG4BlockDisplayAction::~PHG4BlockDisplayAction()
{
  delete m_VisAtt;
}

void PHG4BlockDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
{
  // check if vis attributes exist, if so someone else has set them and we do nothing
  if (m_MyVolume->GetVisAttributes())
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
  m_MyVolume->SetVisAttributes(m_VisAtt);
  return;
}
