#include "PHG4CylinderDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Utils.h>

#include <phparameter/PHParameters.h>

#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

class G4VPhysicalVolume;

using namespace std;

PHG4CylinderDisplayAction::PHG4CylinderDisplayAction(const std::string &name, PHParameters *pars)
  : PHG4DisplayAction(name)
  , m_Params(pars)
  , m_MyVolume(nullptr)
  , m_VisAtt(nullptr)
{
}

PHG4CylinderDisplayAction::~PHG4CylinderDisplayAction()
{
  delete m_VisAtt;
}

void PHG4CylinderDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
{
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
