#include "PHG4CylinderDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <Geant4/G4VPhysicalVolume.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4VisAttributes.hh>

#include <iostream>

using namespace std;
PHG4CylinderDisplayAction::PHG4CylinderDisplayAction(const std::string &name, PHParameters *pars):
  PHG4DisplayAction(name),
  m_Params(pars)
{}

void PHG4CylinderDisplayAction::ApplyDisplayAction(G4VPhysicalVolume* physvol)
{
  cout << "applying display action for " << GetName() << " physvol: " 
       <<  physvol->GetName() << endl;
  FindVolumes(physvol);

}

bool PHG4CylinderDisplayAction::CheckVolume(G4VPhysicalVolume *physvol)
{
  if (physvol == m_MyVolume)
  {
    return true;
  }
  return false;
}

void PHG4CylinderDisplayAction::ApplyVisAttributes(G4VPhysicalVolume *physvol)
{
  G4LogicalVolume* myvol = physvol->GetLogicalVolume();
  G4VisAttributes* visatt = const_cast<G4VisAttributes*> (myvol->GetVisAttributes());
  if (visatt) 
  {
    cout << "visibility: " << visatt->IsVisible() 
	 << " Red: " << visatt->GetColor().GetRed() << endl;
    visatt->SetColor(G4Colour::Magenta());
  }
  else
  {
    visatt = new G4VisAttributes();
    visatt->SetVisibility(true);
    visatt->SetColor(G4Colour::Red());
    myvol->SetVisAttributes(visatt);
  }
  return;
}

