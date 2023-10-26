#ifndef QA_QAG4UTIL_H
#define QA_QAG4UTIL_H

/*!
 * \file QAG4Util.h
 * \brief some common utility functions used in detector specific evaluation modules
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <g4main/PHG4Hit.h>

#include <cmath>

// utility functions
namespace QAG4Util
{
  /// square
  template <class T>
  inline constexpr T square(T x)
  {
    return x * x;
  }

  /// radius
  template <class T>
  inline constexpr T get_r(T x, T y)
  {
    return std::sqrt(square(x) + square(y));
  }

  /// angle difference between [-pi, pi[
  template <class T>
  inline const T delta_phi(T phi1, T phi2)
  {
    auto out = phi1 - phi2;
    while (out >= M_PI) out -= 2 * M_PI;
    while (out < -M_PI) out += 2 * M_PI;
    return out;
  }

  /// calculate the best average of member function called on all members in collection, extrapolated at a given radius
  /** each hit is weighted by its deposited energy */
  template <float (PHG4Hit::*accessor)(int) const>
  float interpolate(std::set<PHG4Hit*> hits, float rextrap)
  {
    // calculate all terms needed for the interpolation
    // need to use double everywhere here due to numerical divergences
    double sw = 0;
    double swr = 0;
    double swr2 = 0;
    double swx = 0;
    double swrx = 0;

    bool valid(false);
    for (const auto& hit : hits)
    {
      const double x0 = (hit->*accessor)(0);
      const double x1 = (hit->*accessor)(1);
      if (std::isnan(x0) || std::isnan(x1)) continue;

      const double w = hit->get_edep();
      if (w < 0) continue;

      valid = true;
      const double r0 = get_r(hit->get_x(0), hit->get_y(0));
      const double r1 = get_r(hit->get_x(1), hit->get_y(1));

      sw += w * 2;
      swr += w * (r0 + r1);
      swr2 += w * (square(r0) + square(r1));
      swx += w * (x0 + x1);
      swrx += w * (r0 * x0 + r1 * x1);
    }

    if (!valid) return NAN;

    const auto alpha = (sw * swrx - swr * swx);
    const auto beta = (swr2 * swx - swr * swrx);
    const auto denom = (sw * swr2 - square(swr));

    return (alpha * rextrap + beta) / denom;
  }

}  // namespace QAG4Util

#endif
