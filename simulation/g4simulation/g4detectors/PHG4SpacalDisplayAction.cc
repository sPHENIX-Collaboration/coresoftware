#include "PHG4SpacalDisplayAction.h"

#include "PHG4CylinderGeom_Spacalv3.h"

#include <g4main/PHG4Utils.h>

#include <TSystem.h>

#include <Geant4/G4Color.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VPhysicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>

#include <iostream>

using namespace std;

PHG4SpacalDisplayAction::PHG4SpacalDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
  , m_Geom(nullptr)
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

void PHG4SpacalDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
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
      visatt->SetVisibility( m_Geom->is_azimuthal_seg_visible() or m_Geom->is_virualize_fiber());
      visatt->SetForceSolid(not m_Geom->is_virualize_fiber());
    }
    else if (it.second == "Divider")
    {
      visatt->SetColor(.8, 1, .8, .3);
      visatt->SetVisibility( m_Geom->is_azimuthal_seg_visible());
    }
    else if (it.second == "Fiber")
    {
      PHG4Utils::SetColour(visatt, "G4_POLYSTYRENE");
      visatt->SetVisibility( m_Geom->is_virualize_fiber());
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
      visatt->SetVisibility( m_Geom->is_azimuthal_seg_visible() or m_Geom->is_virualize_fiber());
      visatt->SetForceSolid(not m_Geom->is_virualize_fiber());
    }
    else if (it.second == "Sector")
    {
      visatt->SetColor(.5, .9, .5, .5);
      visatt->SetVisibility( m_Geom->is_azimuthal_seg_visible() or m_Geom->is_virualize_fiber());
      visatt->SetForceSolid(false);
      visatt->SetForceWireframe(true);
    }
    else if (it.second == "SpacalCylinder")
    {
      PHG4Utils::SetColour(visatt, "W_Epoxy");
      visatt->SetForceSolid((not m_Geom->is_virualize_fiber()) and (not m_Geom->is_azimuthal_seg_visible()));
    }
    else if (it.second == "Wall")
    {
      visatt->SetColor(.5, .9, .5, .1);
      visatt->SetVisibility(m_Geom->is_azimuthal_seg_visible());
    }
    else
    {
      cout << "did not assing color to " << it.first->GetName()
           << " under " << it.second << endl;
      gSystem->Exit(1);
    }
    logvol->SetVisAttributes(visatt);
  }
  return;
}
