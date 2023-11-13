// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EVENTPLANEINFO_H
#define EVENTPLANEINFO_H

#include <phool/PHObject.h>

#include <cmath>
#include <iostream>

class Eventplaneinfo : public PHObject
{
 public:
  ~Eventplaneinfo() override {}

  void identify(std::ostream& os = std::cout) const override
  {
    os << "Eventplaneinfo base class" << std::endl;
  }

  PHObject* CloneMe() const override { return nullptr; }

  virtual void set_qvector(std::vector<std::pair<double, double>> /*Qvec*/) { return; }
  virtual std::pair<double, double> get_qvector(int /*order*/) const { return std::make_pair(NAN, NAN); }
  virtual double get_psi(int /*order*/) const { return NAN; }

 protected:
  Eventplaneinfo() {}

 private:
  ClassDefOverride(Eventplaneinfo, 1);
};

#endif
