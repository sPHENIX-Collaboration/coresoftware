
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
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/detail/VoidPropagatorComponents.hpp>

#include <Acts/Utilities/Result.hpp>

#include <trackbase/ActsGeometry.h>

class SvtxTrack;
class SvtxVertex;
class SvtxVertexMap;
class SvtxTrackState;

class ActsPropagator
{
 public:
  using BoundTrackParam = const Acts::BoundTrackParameters;
  /// Return type of std::pair<path length, parameters>
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

  /// Helper functions for creating needed input for track propagation
  /// functions below
  SurfacePtr makeVertexSurface(const SvtxVertex* vertex);
  SurfacePtr makeVertexSurface(const Acts::Vector3& vertex);
  BoundTrackParam makeTrackParams(SvtxTrack* track, SvtxVertexMap* vertexMap);
  BoundTrackParam makeTrackParams(SvtxTrackState* state, int trackCharge,
                                  SurfacePtr surf);

  /// The return type is an Acts::Result of a std::pair, where the pair is
  /// a path length and the track parameters at the surface in units of mm
  /// and GeV. For an example of how to unpack this, see
  /// PHActsTrackProjection::propagateTrack and
  /// PHActsTrackProjection::updateSvtxTrack
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
  void setConstFieldValue(float field) { m_fieldval = field; }
  void constField() { m_constField = true; }
  void setOverstepLimit(const double overstep) { m_overstepLimit = overstep; }
  SphenixPropagator makePropagator();
  FastPropagator makeFastPropagator();

 private:
  void printTrackParams(const Acts::BoundTrackParameters& params);

  int m_verbosity = 0;

  bool m_constField = false;

  ActsGeometry* m_geometry = nullptr;

  float m_fieldval = 1.4 * Acts::UnitConstants::T;

  /// Default Acts limit
  float m_overstepLimit = 0.01 * Acts::UnitConstants::cm;  // sphenix units cm
};

#endif
