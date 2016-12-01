#include "PHG4SiliconTrackerParameterisation.h"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

#include <boost/format.hpp>

PHG4SiliconTrackerStripParameterisation::PHG4SiliconTrackerStripParameterisation(const unsigned int &ny, const unsigned int &nz, const double &dy, const double &dz)
    : G4VPVParameterisation()
{
  const double offsety = (double)ny * dy * 0.5;
  const double offsetz = (double)nz * dz * 0.5;

  int icopy = 0;
  for (unsigned int iy=0; iy<ny; iy++)
    for (unsigned int iz=0; iz<nz; iz++)
      {
        fXStrip[icopy] = 0.0;
        fYStrip[icopy] = ((double)iy+0.5) * dy - offsety;
        fZStrip[icopy] = ((double)iz+0.5) * dz - offsetz;

        icopy++;
      }
}

PHG4SiliconTrackerStripParameterisation::~PHG4SiliconTrackerStripParameterisation()
{}

void PHG4SiliconTrackerStripParameterisation::ComputeTransformation(const G4int icopy, G4VPhysicalVolume *physVol) const
  {
    physVol->SetTranslation(G4ThreeVector(fXStrip[icopy], fYStrip[icopy], fZStrip[icopy]));
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PHG4SiliconTrackerFPHXParameterisation::PHG4SiliconTrackerFPHXParameterisation(const double &offsetx, const double &offsety, const double &offsetz, const double &dz, const int &ncopy)
    : G4VPVParameterisation()
{
  for (G4int icopy=0; icopy<ncopy; icopy++)
    {
      fXFPHX[icopy] = offsetx;
      fYFPHX[icopy] = offsety;
      fZFPHX[icopy] = offsetz + (double)icopy*dz;
    }
}

PHG4SiliconTrackerFPHXParameterisation::~PHG4SiliconTrackerFPHXParameterisation()
{}

void PHG4SiliconTrackerFPHXParameterisation::ComputeTransformation(const G4int icopy, G4VPhysicalVolume *physVol) const
  {
    physVol->SetTranslation(G4ThreeVector(fXFPHX[icopy], fYFPHX[icopy], fZFPHX[icopy]));
  }
