#ifndef G4MAIN_PHG4IONGUN_H
#define G4MAIN_PHG4IONGUN_H

#include "PHG4PrimaryGeneratorAction.h"

#include <string>

class PHG4IonGun : public PHG4PrimaryGeneratorAction
{
 public:
  PHG4IonGun();
  virtual ~PHG4IonGun() {}

  virtual void GeneratePrimaries(G4Event* anEvent);
  void SetA(const int a) { A = a; }
  void SetZ(const int z) { Z = z; }
  void SetCharge(const int c);
  void ExcitEnergy(const double e) { excitEnergy = e; }
  void SetMom(const double px, const double py, const double pz);
  void Print(const std::string &what = "ALL") const;

 protected:
  int A;
  int Z;
  double mom[3];
  double ioncharge;
  double excitEnergy;
};

#endif
