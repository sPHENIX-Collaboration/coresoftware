/*!
 * \file TrackSeedHelper.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@lanl.gov>
 */

#include "TrackSeedHelper.h"
#include "TrackSeed.h"

#include <trackbase/TrackFitUtils.h>

namespace
{

  //! convenience square method
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }

  std::pair<float, float> findRoot(TrackSeed const* seed)
  {
    /**
    * We need to determine the closest point on the circle to the origin
    * since we can't assume that the track originates from the origin
    * The eqn for the circle is (x-X0)^2+(y-Y0)^2=R^2 and we want to
    * minimize d = sqrt((0-x)^2+(0-y)^2), the distance between the
    * origin and some (currently, unknown) point on the circle x,y.
    *
    * Solving the circle eqn for x and substituting into d gives an eqn for
    * y. Taking the derivative and setting equal to 0 gives the following
    * two solutions. We take the smaller solution as the correct one, as
    * usually one solution is wildly incorrect (e.g. 1000 cm)
    */
    const auto qOverR = seed->get_qOverR();
    const auto X0 = seed->get_X0();
    const auto Y0 = seed->get_Y0();

    const float R = std::abs(1./qOverR);
    const double miny = (std::sqrt(square(X0) * square(R) * square(Y0) + square(R) * pow(Y0, 4)) + square(X0) * Y0 + pow(Y0, 3)) / (square(X0) + square(Y0));
    const double miny2 = (-std::sqrt(square(X0) * square(R) * square(Y0) + square(R) * pow(Y0, 4)) + square(X0) * Y0 + pow(Y0, 3)) / (square(X0) + square(Y0));

    const double minx = std::sqrt(square(R) - square(miny - Y0)) + X0;
    const double minx2 = -std::sqrt(square(R) - square(miny2 - Y0)) + X0;

    /// Figure out which of the two roots is actually closer to the origin
    const float x = (std::abs(minx) < std::abs(minx2)) ? minx : minx2;
    const float y = (std::abs(miny) < std::abs(miny2)) ? miny : miny2;
    return {x, y};
  }

}

//____________________________________________________________________________________
float TrackSeedHelper::get_phi(TrackSeed const* seed, const TrackSeedHelper::position_map_t& positions)
{
  const auto X0 = seed->get_X0();
  const auto Y0 = seed->get_Y0();
  const auto [x, y] = findRoot(seed);

  // This is the angle of the tangent to the circle
  // The argument is the slope of the tangent (inverse of slope of radial line at tangent)
  float phi = std::atan2(-1 * (X0 - x), (Y0 - y));
  Acts::Vector3 pos0 = positions.find(*(seed->begin_cluster_keys()))->second;
  Acts::Vector3 pos1 = positions.find(*std::next(seed->begin_cluster_keys()))->second;

  // we need to know if the track proceeds clockwise or CCW around the circle
  double dx0 = pos0(0) - X0;
  double dy0 = pos0(1) - Y0;
  double phi0 = std::atan2(dy0, dx0);
  double dx1 = pos1(0) - X0;
  double dy1 = pos1(1) - Y0;
  double phi1 = std::atan2(dy1, dx1);
  double dphi = phi1 - phi0;

  // need to deal with the switch from -pi to +pi at phi = 180 degrees
  // final phi - initial phi must be < 180 degrees for it to be a valid track
  if (dphi > M_PI)
  {
    dphi -= 2.0 * M_PI;
  }
  if (dphi < -M_PI)
  {
    dphi += M_PI;
  }

  // whether we add 180 degrees depends on the angle of the bend
  if (dphi < 0)
  {
    phi += M_PI;
    if (phi > M_PI)
    {
      phi -= 2. * M_PI;
    }
  }

  return phi;
}

//____________________________________________________________________________________
float TrackSeedHelper::get_phi_fastsim(TrackSeed const* seed)
{
  const auto [x, y] = findRoot(seed);
  return std::atan2(-(seed->get_X0()-x), (seed->get_Y0() - y));
}

//____________________________________________________________________________________
void TrackSeedHelper::circleFitByTaubin(
  TrackSeed* seed,
  const TrackSeedHelper::position_map_t& positions,
  uint8_t startLayer,
  uint8_t endLayer)
{
  TrackFitUtils::position_vector_t positions_2d;
  for( auto key_iter = seed->begin_cluster_keys(); key_iter != seed->end_cluster_keys(); ++key_iter )
  {
    const auto& key(*key_iter);
    const auto layer = TrkrDefs::getLayer(key);
    if (layer < startLayer or layer > endLayer)
    {
      continue;
    }

    const auto iter = positions.find(key);

    /// you supplied the wrong key...
    if (iter == positions.end())
    {
      continue;
    }

    // add to 2d position list
    const Acts::Vector3& pos = iter->second;
    positions_2d.emplace_back(pos.x(), pos.y());
  }

  // do the fit
  const auto [r, x0, y0] = TrackFitUtils::circle_fit_by_taubin(positions_2d);
  float qOverR = 1./r;

  /// Set the charge
  const auto& firstpos = positions_2d.at(0);
  const auto& secondpos = positions_2d.at(1);

  const auto firstphi = atan2(firstpos.second, firstpos.first);
  const auto secondphi = atan2(secondpos.second, secondpos.first);
  auto dphi = secondphi - firstphi;
  if (dphi > M_PI)
  {
    dphi = 2.*M_PI-dphi;
  }

  if (dphi < -M_PI)
  {
    dphi = 2*M_PI + dphi;
  }
  if (dphi > 0)
  {
    qOverR *= -1;
  }

  // assign
  seed->set_X0(x0);
  seed->set_Y0(y0);
  seed->set_qOverR(qOverR);
}

//____________________________________________________________________________________
void TrackSeedHelper::lineFit(
  TrackSeed* seed,
  const TrackSeedHelper::position_map_t& positions,
  uint8_t startLayer,
  uint8_t endLayer)
{
  TrackFitUtils::position_vector_t positions_2d;
  for( auto key_iter = seed->begin_cluster_keys(); key_iter != seed->end_cluster_keys(); ++key_iter )
  {
    const auto& key(*key_iter);
    const auto layer = TrkrDefs::getLayer(key);
    if (layer < startLayer or layer > endLayer)
    {
      continue;
    }

    const auto iter = positions.find(key);

    /// The wrong key was supplied...
    if (iter == positions.end())
    {
      continue;
    }

    // store (r,z)
    const Acts::Vector3& pos = iter->second;
    positions_2d.emplace_back(std::sqrt(square(pos.x()) + square(pos.y())), pos.z());
  }

  // do the fit
  const auto [slope, intercept] = TrackFitUtils::line_fit(positions_2d);

  // assign
  seed->set_slope(slope);
  seed->set_Z0(intercept);
}

//____________________________________________________________________________________
float TrackSeedHelper::get_x(TrackSeed const* seed)
{
  return findRoot(seed).first;
}

//____________________________________________________________________________________
float TrackSeedHelper::get_y(TrackSeed const* seed)
{
  return findRoot(seed).second;
}

//____________________________________________________________________________________
float TrackSeedHelper::get_z(TrackSeed const* seed)
{
  return seed->get_Z0();
}

//____________________________________________________________________________________
Acts::Vector3 TrackSeedHelper::get_xyz(TrackSeed const* seed)
{
  const auto [x,y] = findRoot(seed);
  return {x,y,seed->get_Z0()};
}
