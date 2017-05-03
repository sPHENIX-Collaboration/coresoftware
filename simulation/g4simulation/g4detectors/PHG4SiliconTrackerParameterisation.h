#ifndef PHG4SiliconTrackerParameterisation_H
#define PHG4SiliconTrackerParameterisation_H 1

#include <Geant4/G4VPVParameterisation.hh>
#include <Geant4/globals.hh>

class G4VPhysicalVolume;

/*
 * Strip location
 */
class PHG4SiliconTrackerStripParameterisation : public G4VPVParameterisation
{
 public:
  PHG4SiliconTrackerStripParameterisation(const unsigned int ny, const unsigned int nz, const double dy, const double dz);
  virtual ~PHG4SiliconTrackerStripParameterisation() {}
  virtual void ComputeTransformation(const G4int icopy, G4VPhysicalVolume *physVol) const;

 private:
  G4double fXStrip[128 * 2 * 8];
  G4double fYStrip[128 * 2 * 8];
  G4double fZStrip[128 * 2 * 8];
};

/*
 * FPHX location
 */
class PHG4SiliconTrackerFPHXParameterisation : public G4VPVParameterisation
{
 public:
  PHG4SiliconTrackerFPHXParameterisation(const double offsetx, const double offsety, const double offsetz, const double dz, const int ncopy);
  virtual ~PHG4SiliconTrackerFPHXParameterisation() {}
  virtual void ComputeTransformation(const G4int icopy, G4VPhysicalVolume *physVol) const;

 private:
  G4double fXFPHX[20];
  G4double fYFPHX[20];
  G4double fZFPHX[20];
};

#endif
