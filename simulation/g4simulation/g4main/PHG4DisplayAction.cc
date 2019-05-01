#include "PHG4DisplayAction.h"

#include <Geant4/G4VPhysicalVolume.hh>
#include <Geant4/G4LogicalVolume.hh>

using namespace std;

void PHG4DisplayAction::FindVolumes(G4VPhysicalVolume *physvol)
{
  if (CheckVolume(physvol))
  {
    ApplyVisAttributes(physvol);
  }
  G4LogicalVolume* logvol = physvol->GetLogicalVolume();
  int nDaughters = logvol->GetNoDaughters();
  if (nDaughters > 0)
  {
    for (int i = 0; i<nDaughters; ++i)
    {
      G4VPhysicalVolume* daughtervol = logvol->GetDaughter(i);
      //G4LogicalVolume* logicdaughter = daughtervol->GetLogicalVolume();
    
      FindVolumes(daughtervol);
    }
  }
}
