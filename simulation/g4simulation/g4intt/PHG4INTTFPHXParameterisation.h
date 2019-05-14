// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4INTT_PHG4INTTPARAMETERISATION_H
#define G4INTT_PHG4INTTPARAMETERISATION_H

#include <Geant4/G4VPVParameterisation.hh>

class G4VPhysicalVolume;

/*
 * FPHX location
 */
class PHG4INTTFPHXParameterisation : public G4VPVParameterisation
{
 public:
  PHG4INTTFPHXParameterisation(const double offsetx, const double offsety, const double offsetz, const double dz, const int ncopy);
  virtual ~PHG4INTTFPHXParameterisation() {}
  virtual void ComputeTransformation(const G4int icopy, G4VPhysicalVolume *physVol) const;

 private:
  G4double fXFPHX[20];
  G4double fYFPHX[20];
  G4double fZFPHX[20];
};

#endif
