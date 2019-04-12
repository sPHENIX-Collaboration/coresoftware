#ifndef PHG4PARTICLEV3_H__
#define PHG4PARTICLEV3_H__

#include "PHG4Particlev2.h"

class PHG4Particlev3: public PHG4Particlev2
{
 public:
  PHG4Particlev3();
//  PHG4Particlev3(const std::string &name, const int pid, const double px, const double py, const double pz);
  PHG4Particlev3(const PHG4Particle *in);

  virtual ~PHG4Particlev3() {}
  bool isIon() const {return true;}
  void set_A(const int a) { A = a; }
  int  get_A() const { return A; }
  void set_Z(const int z) { Z = z; }
  int get_Z() const { return Z; }
  void set_NumCharge(const int c); // number of charges - gets converted to charge
  void set_IonCharge(const double ch) {ioncharge=ch;}
  double get_IonCharge() const {return ioncharge;}
  void set_ExcitEnergy(const double e) { excitEnergy = e; }
  double get_ExcitEnergy() const { return excitEnergy; }
  void identify(std::ostream& os = std::cout) const;

 protected:
  int A;
  int Z;
  double ioncharge;
  double excitEnergy;

  ClassDef(PHG4Particlev3,1)
};

#endif
