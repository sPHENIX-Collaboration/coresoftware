// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4IONGUN_H
#define G4MAIN_PHG4IONGUN_H

#include "PHG4ParticleGeneratorBase.h"

#include <string>

class PHG4Particle;
class PHCompositeNode;

class PHG4IonGun : public PHG4ParticleGeneratorBase
{
 public:
  PHG4IonGun(const std::string &name = "PHG4IONGUN");
  virtual ~PHG4IonGun() {}
  int process_event(PHCompositeNode *topNode);
  void SetA(const int a) { A = a; }
  void SetZ(const int z) { Z = z; }
  void SetCharge(const int c);
  void ExcitEnergy(const double e) { excitEnergy = e; }
  void SetMom(const double px, const double py, const double pz);
  void Print(const std::string &what = "ALL") const;

 private:
  void UpdateParticle();
  PHG4Particle *ion;
  int A;
  int Z;
  double mom[3];
  int ioncharge;
  double excitEnergy;
};

#endif
