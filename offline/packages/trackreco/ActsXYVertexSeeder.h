// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TRACKRECO_ACTSXYVERTEXSEEDER_H
#define TRACKRECO_ACTSXYVERTEXSEEDER_H

/*!
 *  \file   ActsXYVertexSeeder.h
 *  \brief  Vertex seed finder for the Acts IterativeVertexFinder that
 *          estimates the transverse seed position from the tracks.
 *
 *  The Acts seed finders (ZScan, TrackDensity) only determine z and place
 *  the seed at the beam spot x,y. The IVF then keeps only tracks within
 *  significanceCutSeeding * sqrt(sigma_d0^2 + sigma_z0^2) of the seed, which
 *  for precise silicon tracks is a window of ~100 um. Collisions displaced
 *  from the beam spot by more than that lose tracks or are split.
 *
 *  This seeder takes the z estimate of an inner seed finder (normally
 *  ZScan), selects the tracks in a z window around it and fits their common
 *  point in x,y (straight-line approximation of each track around its PCA,
 *  weighted by sigma_d0, with the beam spot as a weak prior), then z, with
 *  iterative outlier rejection. If fewer than two tracks survive, the inner
 *  seed is returned unchanged.
 */

#include <Acts/Utilities/Logger.hpp>
#include <Acts/Vertexing/IVertexFinder.hpp>
#include <Acts/Vertexing/TrackAtVertex.hpp>

#include <memory>

class ActsXYVertexSeeder final : public Acts::IVertexFinder
{
 public:
  struct Config
  {
    /// seed finder providing the z estimate. Its x,y are ignored
    std::shared_ptr<const Acts::IVertexFinder> zSeeder;

    /// extracts the track parameters from an InputTrack
    Acts::InputTrack::Extractor extractParameters;

    /// half width of the z window around the current seed z used to select
    /// tracks (Acts units, mm). Must cover the z shift of the tracks' PCA
    /// to the beam axis, ~ cot(theta) * transverse displacement
    double zWindow = 2.;

    /// tracks with (d0/sigma_d0)^2 + (dz/sigma_z0)^2 above this w.r.t. the
    /// fitted seed are dropped from the seed estimate (one at a time)
    double maxTrackChi2 = 36.;

    /// number of select-and-fit passes
    int nIterations = 3;
  };

  /// state for fulfilling the interface
  struct State
  {
  };

  explicit ActsXYVertexSeeder(Config cfg,
                              std::unique_ptr<const Acts::Logger> logger =
                                  Acts::getDefaultLogger("ActsXYVertexSeeder", Acts::Logging::INFO));

  Acts::Result<std::vector<Acts::Vertex>> find(
      const std::vector<Acts::InputTrack>& trackVector,
      const Acts::VertexingOptions& vertexingOptions,
      Acts::IVertexFinder::State& state) const override;

  Acts::IVertexFinder::State makeState(
      const Acts::MagneticFieldContext& /*mctx*/) const override
  {
    return Acts::IVertexFinder::State{State{}};
  }

  void setTracksToRemove(
      Acts::IVertexFinder::State& /*state*/,
      const std::vector<Acts::InputTrack>& /*removedTracks*/) const override
  {
  }

 private:
  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
  const Acts::Logger& logger() const { return *m_logger; }
};

#endif  // TRACKRECO_ACTSXYVERTEXSEEDER_H
