// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EVENTPLANEINFOV2_H
#define EVENTPLANEINFOV2_H

#include "Eventplaneinfo.h"

#include <cstddef>  // for size_t
#include <iostream>
#include <limits>
#include <utility>  // for pair, make_pair
#include <vector>

class PHObject;

class Eventplaneinfov2 : public Eventplaneinfo
{
 public:
  Eventplaneinfov2() = default;
  ~Eventplaneinfov2() override = default;

  Eventplaneinfov2(const Eventplaneinfov2&) = default;
  Eventplaneinfov2& operator=(const Eventplaneinfov2&) = default;
  Eventplaneinfov2(Eventplaneinfov2&&) = default;
  Eventplaneinfov2& operator=(Eventplaneinfov2&&) = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = Eventplaneinfov2(); }
  PHObject* CloneMe() const override { return new Eventplaneinfov2(*this); }

  void set_qvector(std::vector<std::pair<double, double>> Qvec) override { mQvec = Qvec; }
  void set_qvector_raw(const std::vector<std::pair<double, double>>& Qvec) override { mQvec_raw = Qvec; }
  void set_qvector_recentered(const std::vector<std::pair<double, double>>& Qvec) override { mQvec_recentered = Qvec; }
  void set_shifted_psi(std::vector<double> Psi_Shifted) override { mPsi_Shifted = Psi_Shifted; }
  std::pair<double, double> get_qvector(int order) const override { return std::make_pair(mQvec[order - 1].first, mQvec[order - 1].second); }
  std::pair<double, double> get_qvector_raw(int order) const override { return std::make_pair(mQvec_raw[order - 1].first, mQvec_raw[order - 1].second); }
  std::pair<double, double> get_qvector_recentered(int order) const override { return std::make_pair(mQvec_recentered[order - 1].first, mQvec_recentered[order - 1].second); }
  void set_ring_qvector(std::vector<std::vector<std::pair<double, double>>> Qvec) override { ring_Qvec = Qvec; }
  std::pair<double, double> get_ring_qvector(int ring_index, int order) const override { return ring_Qvec[ring_index][order - 1]; }
  double get_ring_psi(int ring_index, int order) const override {return GetPsi(ring_Qvec[ring_index][order - 1].first,ring_Qvec[ring_index][order - 1].second,order);}
  double GetPsi(double Qx, double Qy, unsigned int order) const override;
  double get_psi(int order) const override { return GetPsi(mQvec[order - 1].first, mQvec[order - 1].second, order);}
  double get_shifted_psi(int order) const override { return mPsi_Shifted[order - 1]; }
     
 private:
  std::vector<std::pair<double, double>> mQvec;
  std::vector<std::pair<double, double>> mQvec_raw;
  std::vector<std::pair<double, double>> mQvec_recentered;
  std::vector<double> mPsi_Shifted;
  std::vector<std::vector<std::pair<double, double>>> ring_Qvec;
  ClassDefOverride(Eventplaneinfov2, 1);
};

#endif
     
