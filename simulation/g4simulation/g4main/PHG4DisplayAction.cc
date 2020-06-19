#include "PHG4DisplayAction.h"

#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4VPhysicalVolume.hh>

int PHG4DisplayAction::FindVolumes(G4VPhysicalVolume* physvol)
{
  int iret = 0;
  if (int chk = CheckVolume(physvol))
  {
    ApplyVisAttributes(physvol);
    if (chk == CheckReturnCodes::ABORT)
    {
      iret = -1;
      return iret;
    }
  }
  G4LogicalVolume* logvol = physvol->GetLogicalVolume();
  int nDaughters = logvol->GetNoDaughters();
  if (nDaughters > 0)
  {
    for (int i = 0; i < nDaughters; ++i)
    {
      G4VPhysicalVolume* daughtervol = logvol->GetDaughter(i);
      //G4LogicalVolume* logicdaughter = daughtervol->GetLogicalVolume();

      if (FindVolumes(daughtervol))
      {
        iret = -1;
        return iret;
      }
    }
  }
  return iret;
}
