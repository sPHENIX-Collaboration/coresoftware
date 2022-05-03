// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4INTT_PHG4INTTPARAMETERISATION_H
#define G4INTT_PHG4INTTPARAMETERISATION_H

#include <Geant4/G4Types.hh>  // for G4double, G4int
#include <Geant4/G4VPVParameterisation.hh>

class G4VPhysicalVolume;

/*
 * FPHX location
 */
class PHG4InttFPHXParameterisation : public G4VPVParameterisation
{
 public:
  PHG4InttFPHXParameterisation(const double offsetx, const double offsety, const double offsetz, const double dz, const int ncopy);
  ~PHG4InttFPHXParameterisation() override {}
  void ComputeTransformation(const G4int icopy, G4VPhysicalVolume *physVol) const override;

 private:
  G4double fXFPHX[20];
  G4double fYFPHX[20];
  G4double fZFPHX[20];
};

#endif
