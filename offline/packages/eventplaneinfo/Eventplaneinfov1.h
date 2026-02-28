#ifndef EVENTPLANEINFO_EVENTPLANEINFOV1_H
#define EVENTPLANEINFO_EVENTPLANEINFOV1_H

#include "Eventplaneinfo.h"

#include <cstddef>  // for size_t
#include <iostream>
#include <limits>
#include <utility>  // for pair, make_pair
#include <vector>

class PHObject;

class Eventplaneinfov1 : public Eventplaneinfo
{
 public:
  Eventplaneinfov1() = default;
  ~Eventplaneinfov1() override = default;

  Eventplaneinfov1(const Eventplaneinfov1&) = default;
  Eventplaneinfov1& operator=(const Eventplaneinfov1&) = default;
  Eventplaneinfov1(Eventplaneinfov1&&) = default;
  Eventplaneinfov1& operator=(Eventplaneinfov1&&) = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = Eventplaneinfov1(); }
  PHObject* CloneMe() const override { return new Eventplaneinfov1(*this); }

  void set_qvector(const std::vector<std::pair<double, double>>& Qvec) override { mQvec = Qvec; }
  void set_qvector_raw(const std::vector<std::pair<double, double>>& Qvec) override { mQvec_raw = Qvec; }
  void set_qvector_recentered(const std::vector<std::pair<double, double>>& Qvec) override { mQvec_recentered = Qvec; }
  void set_shifted_psi(const std::vector<double>& Psi_Shifted) override { mPsi_Shifted = Psi_Shifted; }
  std::pair<double, double> get_qvector(int order) const override { return safe_qvec(mQvec, order); }
  std::pair<double, double> get_qvector_raw(int order) const override { return safe_qvec(mQvec_raw, order); }
  std::pair<double, double> get_qvector_recentered(int order) const override { return safe_qvec(mQvec_recentered, order); }
  void set_ring_qvector(const std::vector<std::vector<std::pair<double, double>>>& Qvec) override { ring_Qvec = Qvec; }
  std::pair<double, double> get_ring_qvector(int ring_index, int order) const override
  {
    if (ring_index < 0 || static_cast<size_t>(ring_index) >= ring_Qvec.size())
    {
      return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }
    return safe_qvec(ring_Qvec[ring_index], order);
  }
  double get_ring_psi(int ring_index, int order) const override
  {
    auto q = get_ring_qvector(ring_index, order);
    return GetPsi(q.first, q.second, static_cast<unsigned int>(order));
  }

  double GetPsi(double Qx, double Qy, unsigned int order) const override;
  double get_psi(int order) const override
  {
    auto q = get_qvector(order);
    return GetPsi(q.first, q.second, static_cast<unsigned int>(order));
  }
  double get_shifted_psi(int order) const override
  {
    if (order <= 0 || static_cast<size_t>(order) > mPsi_Shifted.size())
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return mPsi_Shifted[order - 1];
  }

 private:
  static std::pair<double, double> safe_qvec(const std::vector<std::pair<double, double>>& v, int order)
  {
    if (order <= 0 || static_cast<size_t>(order) > v.size())
    {
      return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }
    return v[order - 1];
  }

  std::vector<std::pair<double, double>> mQvec;
  std::vector<std::pair<double, double>> mQvec_raw;
  std::vector<std::pair<double, double>> mQvec_recentered;
  std::vector<double> mPsi_Shifted;
  std::vector<std::vector<std::pair<double, double>>> ring_Qvec;
  ClassDefOverride(Eventplaneinfov1, 1);
};

#endif
