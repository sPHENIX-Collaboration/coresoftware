// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PARTICLEV1_H
#define G4MAIN_PHG4PARTICLEV1_H

#include "PHG4Particle.h"

#include <iostream>
#include <string>

class PHG4Particlev1 : public PHG4Particle
{
 public:
  PHG4Particlev1();
  PHG4Particlev1(const std::string &name, const int pid, const double px, const double py, const double pz);
  PHG4Particlev1(const PHG4Particle *in);

  ~PHG4Particlev1() override {}

  void identify(std::ostream &os = std::cout) const override;

  int get_pid() const override { return fpid; }
  std::string get_name() const override { return fname; }
  double get_px() const override { return fpx; }
  double get_py() const override { return fpy; }
  double get_pz() const override { return fpz; }

  int get_barcode() const override { return barcode; }

  void set_name(const std::string &name) override { fname = name; }
  void set_pid(const int i) override { fpid = i; }
  void set_px(const double x) override { fpx = x; }
  void set_py(const double x) override { fpy = x; }
  void set_pz(const double x) override { fpz = x; }

  void set_barcode(const int bcd) override { barcode = bcd; }


 protected:
  std::string fname;
  int fpid;
  double fpx, fpy, fpz;
  int barcode;

  ClassDefOverride(PHG4Particlev1, 1)
};

#endif
