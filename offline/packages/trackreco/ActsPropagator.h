
#ifndef ACTSPROPAGATOR_H
#define ACTSPROPAGATOR_H

#include <trackbase/ActsGeometry.h>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/EventData/TrackParameters.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Acts/Propagator/EigenStepper.hpp>
#pragma GCC diagnostic pop
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Acts/Propagator/Propagator.hpp>
#pragma GCC diagnostic pop
#include <Acts/Propagator/detail/VoidPropagatorComponents.hpp>
#include <Acts/Propagator/Navigator.hpp>

#include <Acts/Utilities/Result.hpp>

#include <trackbase/ActsGeometry.h>

class ActsPropagator
{
 public:
  using BoundTrackParam = const Acts::BoundTrackParameters;
  using BoundTrackParamPair = std::pair<float, BoundTrackParam>;
  using BoundTrackParamResult = Acts::Result<BoundTrackParamPair>;
  using SurfacePtr = std::shared_ptr<const Acts::Surface>;
  using Stepper = Acts::EigenStepper<>;
  using FastPropagator = Acts::Propagator<Stepper>;
  using SphenixPropagator = Acts::Propagator<Stepper, Acts::Navigator>;

  ActsPropagator() {}
  ActsPropagator(ActsGeometry* geometry)
    : m_geometry(geometry)
  {
  }
  ~ActsPropagator() {}

  BoundTrackParamResult propagateTrack(const Acts::BoundTrackParameters& params,
                                       const unsigned int sphenixLayer);

  BoundTrackParamResult propagateTrack(const Acts::BoundTrackParameters& params,
                                       const SurfacePtr& surface);

  /// The following function takes the track parameters at the vertex and
  /// propagates them in isolation to the requested surface, i.e. it does
  /// NOT stop at each layer in the sPHENIX detector on the way to the
  /// target surface
  BoundTrackParamResult propagateTrackFast(const Acts::BoundTrackParameters& params,
                                           const SurfacePtr& surface);

  bool checkLayer(const unsigned int& sphenixlayer,
                  unsigned int& actsvolume,
                  unsigned int& actslayer);
  void verbosity(int verb) { m_verbosity = verb; }

 private:
  SphenixPropagator makePropagator();
  FastPropagator makeFastPropagator();
  void printTrackParams(const Acts::BoundTrackParameters& params);

  int m_verbosity = 0;

  bool m_constField = false;

  ActsGeometry* m_geometry = nullptr;
};

#endif
