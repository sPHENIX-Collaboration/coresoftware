// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EVENTPLANEINFO_EVENTPLANEINFO_H
#define EVENTPLANEINFO_EVENTPLANEINFO_H

#include <phool/PHObject.h>

#include <iostream>
#include <limits>
#include <vector>
#include <utility>

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
  virtual void set_qvector_raw(const std::vector<std::pair<double, double>>& /*Qvec*/) { return; }
  virtual void set_qvector_recentered(const std::vector<std::pair<double, double>>& /*Qvec*/) { return; }
  virtual void set_shifted_psi(std::vector<double> /*Psi_Shifted*/) { return; }
  virtual std::pair<double, double> get_qvector(int /*order*/) const { return std::make_pair(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()); }
  virtual std::pair<double, double> get_qvector_raw(int /*order*/) const { return std::make_pair(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()); }
  virtual std::pair<double, double> get_qvector_recentered(int /*order*/) const { return std::make_pair(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()); }
  virtual double get_psi(int /*order*/) const { return std::numeric_limits<double>::quiet_NaN(); }
  virtual double get_shifted_psi(int /*order*/) const { return std::numeric_limits<double>::quiet_NaN(); }
  virtual double GetPsi(const double /*Qx*/, const double /*Qy*/, const unsigned int /*order*/) const { return std::numeric_limits<double>::quiet_NaN(); }
  virtual void set_ring_qvector(std::vector<std::vector<std::pair<double, double>>> /*RingQvecs*/) { return; }
  virtual std::pair<double, double> get_ring_qvector(int /*rbin*/, int /*order*/) const { return std::make_pair(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()); }
  virtual double get_ring_psi(int /*rbin*/, int /*order*/) const { return std::numeric_limits<double>::quiet_NaN(); }

 protected:
  Eventplaneinfo() = default;

 private:
  ClassDefOverride(Eventplaneinfo, 1);
};

#endif
