#include "PHG4SiliconTrackerParameterisation.h"

#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4VPhysicalVolume.hh>

PHG4SiliconTrackerStripParameterisation::PHG4SiliconTrackerStripParameterisation(const unsigned int ny, const unsigned int nz, const double dy, const double dz)
{
  const double offsety = ny * dy * 0.5;
  const double offsetz = nz * dz * 0.5;
 
  int icopy = 0;
  for (unsigned int iy = 0; iy < ny; iy++)
    for (unsigned int iz = 0; iz < nz; iz++)
    {
      fXStrip[icopy] = 0.0;
      fYStrip[icopy] = (iy + 0.5) * dy - offsety;
      fZStrip[icopy] = (iz + 0.5) * dz - offsetz;

      /*
	std::cout << "      icopy " << icopy
		  << " iy " << iy << " iz " << iz
		  << " offsety " << offsety
		  << " offsetz " << offsetz
		  << " fXStrip[icopy] " << fXStrip[icopy]
		  << " fYStrip[icopy] " << fYStrip[icopy]
		  << " fZStrip[icopy] " << fZStrip[icopy]
		  << std::endl;
      */

      icopy++;
    }
}

void PHG4SiliconTrackerStripParameterisation::ComputeTransformation(const G4int icopy, G4VPhysicalVolume *physVol) const
{
  physVol->SetTranslation(G4ThreeVector(fXStrip[icopy], fYStrip[icopy], fZStrip[icopy]));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PHG4SiliconTrackerFPHXParameterisation::PHG4SiliconTrackerFPHXParameterisation(const double offsetx, const double offsety, const double offsetz, const double dz, const int ncopy)
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

void PHG4SiliconTrackerFPHXParameterisation::ComputeTransformation(const G4int icopy, G4VPhysicalVolume *physVol) const
{
  physVol->SetTranslation(G4ThreeVector(fXFPHX[icopy], fYFPHX[icopy], fZFPHX[icopy]));
}
