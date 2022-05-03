#include "PHG4ConeDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Utils.h>

#include <phparameter/PHParameters.h>

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <string>

class G4VPhysicalVolume;

PHG4ConeDisplayAction::PHG4ConeDisplayAction(const std::string &name, PHParameters *pars)
  : PHG4DisplayAction(name)
  , m_Params(pars)
  , m_MyVolume(nullptr)
  , m_VisAtt(nullptr)
  , m_Colour(nullptr)
{
}

PHG4ConeDisplayAction::~PHG4ConeDisplayAction()
{
  delete m_VisAtt;
  delete m_Colour;
}

void PHG4ConeDisplayAction::ApplyDisplayAction(G4VPhysicalVolume * /*physvol*/)
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
  if (m_Colour)
  {
    m_VisAtt->SetColour(m_Colour->GetRed(),
                        m_Colour->GetGreen(),
                        m_Colour->GetBlue(),
                        m_Colour->GetAlpha());
    m_VisAtt->SetVisibility(true);
    m_VisAtt->SetForceSolid(true);
  }
  // drawing 200 segments per circle makes it look smoother than default
  m_VisAtt->SetForceLineSegmentsPerCircle(200);
  m_MyVolume->SetVisAttributes(m_VisAtt);
  return;
}

void PHG4ConeDisplayAction::SetColor(const double red, const double green, const double blue, const double alpha)
{
  if (std::isfinite(red) && std::isfinite(green) && std::isfinite(blue) && std::isfinite(alpha))
  {
    m_Colour = new G4Colour(red, green, blue, alpha);
  }
  return;
}
