#include "PHG4SpacalDisplayAction.h"

#include "PHG4CylinderGeom_Spacalv1.h"  // for PHG4CylinderGeom_Spacalv1

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Utils.h>

#include <TSystem.h>

#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4VisAttributes.hh>

#include <iostream>
#include <utility>  // for pair

PHG4SpacalDisplayAction::PHG4SpacalDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4SpacalDisplayAction::~PHG4SpacalDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4SpacalDisplayAction::ApplyDisplayAction(G4VPhysicalVolume * /*physvol*/)
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
    visatt->SetForceWireframe(false);

    m_VisAttVec.push_back(visatt);  // for later deletion
    if (it.second == "AzimuthSegment")
    {
      visatt->SetColor(.1, .1, .1, .5);
      visatt->SetVisibility(m_Geom->is_virualize_fiber());
      visatt->SetForceSolid(false);
    }
    else if (it.second == "Block")
    {
      visatt->SetColor(.3, .3, .3, .3);
      visatt->SetVisibility(m_Geom->is_azimuthal_seg_visible() or m_Geom->is_virualize_fiber());
      visatt->SetForceSolid(!m_Geom->is_virualize_fiber());
    }
    else if (it.second == "Divider")
    {
      visatt->SetColor(.8, 1, .8, .3);
      visatt->SetVisibility(m_Geom->is_azimuthal_seg_visible());
    }
    else if (it.second == "Fiber")
    {
      PHG4Utils::SetColour(visatt, "G4_POLYSTYRENE");
      visatt->SetVisibility(m_Geom->is_virualize_fiber());
      visatt->SetForceSolid(m_Geom->is_virualize_fiber());
    }
    else if (it.second == "FiberCore")
    {
      PHG4Utils::SetColour(visatt, "G4_POLYSTYRENE");
      visatt->SetVisibility(false);
      visatt->SetForceSolid(false);
    }
    else if (it.second == "LightGuide")
    {
      PHG4Utils::SetColour(visatt, m_MaterialMap["LightGuide"]);
      visatt->SetColor(.8, 1, .8, .3);
      visatt->SetVisibility(m_Geom->is_azimuthal_seg_visible() or m_Geom->is_virualize_fiber());
      visatt->SetForceSolid(!m_Geom->is_virualize_fiber());
    }
    else if (it.second == "Sector")
    {
      visatt->SetColor(.5, .9, .5, .5);
      visatt->SetVisibility(m_Geom->is_azimuthal_seg_visible() or m_Geom->is_virualize_fiber());
      visatt->SetForceSolid(false);
      visatt->SetForceWireframe(true);
    }
    else if (it.second == "SpacalCylinder")
    {
      PHG4Utils::SetColour(visatt, "W_Epoxy");
      visatt->SetForceSolid((!m_Geom->is_virualize_fiber()) && (!m_Geom->is_azimuthal_seg_visible()));
    }
    else if (it.second == "Wall")
    {
      visatt->SetColor(.5, .9, .5, .1);
      visatt->SetVisibility(m_Geom->is_azimuthal_seg_visible());
    }
    else if (it.second == "WallProj")
    {
      visatt->SetColor(.5, .9, .5, .2);
      visatt->SetVisibility(m_Geom->is_azimuthal_seg_visible() && (!m_Geom->is_virualize_fiber()));
    }
    else
    {
      std::cout << "did not assing color to " << it.first->GetName()
                << " under " << it.second << std::endl;
      gSystem->Exit(1);
    }
    logvol->SetVisAttributes(visatt);
  }
  return;
}
