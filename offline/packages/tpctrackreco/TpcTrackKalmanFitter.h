// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCTRACKKALMANFITTER_H
#define TPCTRACKRECO_TPCTRACKKALMANFITTER_H

#include "TpcTrackFit.h"

#include <array>
#include <utility>
#include <vector>

class TpcTrackKalmanFitter
{
 public:
  enum StateIndex
  {
    X = 0,
    Y = 1,
    Z = 2,
    Phi = 3,
    QOverPt = 4,
    TanLambda = 5,
    StateDim = 6
  };

  static bool fit(const std::vector<TpcTrackPoint> &points,
                  int charge,
                  const TpcKalmanConfig &config,
                  TpcKalmanResult &result,
                  double mass_gev = 0.13957039);

  static TpcTrackVec3 state_position(const std::array<double, StateDim> &state);
  static TpcTrackVec3 state_momentum(const std::array<double, StateDim> &state);
  static TpcTrackVec3 state_tangent(const std::array<double, StateDim> &state);
  static std::array<double, StateDim> propagation_state(const TpcKalmanResult &fit,
                                                        const TpcTrackVec3 &reference_vertex);
  static std::array<double, StateDim> propagate_state(const std::array<double, StateDim> &state,
                                                      double ds_cm,
                                                      const TpcKalmanConfig &config,
                                                      double mass_gev = 0.13957039);
  static std::pair<double, double> dca_to_vertex(const TpcKalmanResult &fit,
                                                 const TpcTrackVec3 &vertex,
                                                 const TpcKalmanConfig *config = nullptr);
};

#endif
