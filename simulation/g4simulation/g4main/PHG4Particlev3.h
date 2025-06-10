// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PARTICLEV3_H
#define G4MAIN_PHG4PARTICLEV3_H

#include "PHG4Particlev2.h"

#include <iostream>
#include <limits>

class PHG4Particle;

class PHG4Particlev3 : public PHG4Particlev2
{
 public:
  PHG4Particlev3() = default;
  //  PHG4Particlev3(const std::string &name, const int pid, const double px, const double py, const double pz);
  PHG4Particlev3(const PHG4Particle* in);

  ~PHG4Particlev3() override = default;

  PHObject* CloneMe() const override { return new PHG4Particlev3(*this); }

  void identify(std::ostream& os = std::cout) const override;

  bool isIon() const override { return true; }
  void set_A(const int a) override { A = a; }
  int get_A() const override { return A; }
  void set_Z(const int z) override { Z = z; }
  int get_Z() const override { return Z; }
  void set_NumCharge(const int c) override;  // number of charges - gets converted to charge
  void set_IonCharge(const double ch) override { ioncharge = ch; }
  double get_IonCharge() const override { return ioncharge; }
  void set_ExcitEnergy(const double e) override { excitEnergy = e; }
  double get_ExcitEnergy() const override { return excitEnergy; }

 protected:
  int A{0};
  int Z{0};
  double ioncharge{std::numeric_limits<double>::quiet_NaN()};
  double excitEnergy{std::numeric_limits<double>::quiet_NaN()};

  ClassDefOverride(PHG4Particlev3, 1)
};

#endif
