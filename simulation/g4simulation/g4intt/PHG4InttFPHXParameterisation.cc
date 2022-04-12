#include "PHG4InttFPHXParameterisation.h"

#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4Types.hh>  // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PHG4InttFPHXParameterisation::PHG4InttFPHXParameterisation(const double offsetx, const double offsety, const double offsetz, const double dz, const int ncopy)
{
  for (G4int icopy = 0; icopy < ncopy; icopy++)
  {
    fXFPHX[icopy] = offsetx;
    fYFPHX[icopy] = offsety;
    fZFPHX[icopy] = offsetz + icopy * dz;
    /*
    std::cout << "      icopy " << icopy
	      << " offsety " << offsety
	      << " offsetz " << offsetz
	      << " dz " << dz
	      << " fXFPHX[icopy] " << fXFPHX[icopy]
	      << " fYFPHX[icopy] " << fYFPHX[icopy]
	      << " fZFPHX[icopy] " << fZFPHX[icopy]
	      << std::endl;
    */
  }
}

void PHG4InttFPHXParameterisation::ComputeTransformation(const G4int icopy, G4VPhysicalVolume *physVol) const
{
  physVol->SetTranslation(G4ThreeVector(fXFPHX[icopy], fYFPHX[icopy], fZFPHX[icopy]));
}
