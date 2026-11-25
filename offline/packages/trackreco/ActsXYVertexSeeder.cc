#include "ActsXYVertexSeeder.h"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Vertexing/Vertex.hpp>
#include <Acts/Vertexing/VertexingOptions.hpp>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace
{
  /// track quantities used by the seed fit, straight-line approximation
  /// around the track position (its PCA to the reference line)
  struct SeedTrack
  {
    Acts::Vector2 pos;   ///< transverse position of the PCA
    Acts::Vector2 dir;   ///< unit transverse direction
    Acts::Vector2 norm;  ///< unit transverse normal, the d0 direction
    double z = 0;        ///< z of the PCA
    double cotTheta = 0; ///< dz per unit transverse path length
    double wD0 = 0;      ///< 1 / sigma_d0^2
    double wZ0 = 0;      ///< 1 / sigma_z0^2
    bool used = false;

    /// z of the track at its closest transverse approach to v
    double zAt(const Acts::Vector2& v) const
    {
      return z + cotTheta * (v - pos).dot(dir);
    }

    /// signed transverse distance of v from the track
    double d0At(const Acts::Vector2& v) const
    {
      return (v - pos).dot(norm);
    }

    double chi2(const Acts::Vector2& v, double zv) const
    {
      const double d0 = d0At(v);
      const double dz = zAt(v) - zv;
      return d0 * d0 * wD0 + dz * dz * wZ0;
    }
  };

  struct SeedFit
  {
    Acts::Vector2 xy = Acts::Vector2::Zero();
    double z = 0;
    Acts::SquareMatrix2 covXY = Acts::SquareMatrix2::Zero();
    double varZ = 0;
    unsigned int nTracks = 0;
  };

  /// weighted fit of the common point of the used tracks, transverse first
  /// (with optional beam spot prior), then z at that transverse position
  bool fitSeed(const std::vector<SeedTrack>& tracks,
               const Acts::Vector2& priorXY,
               const Acts::SquareMatrix2& priorWeight,
               SeedFit& fit)
  {
    Acts::SquareMatrix2 A = priorWeight;
    Acts::Vector2 b = priorWeight * priorXY;
    unsigned int n = 0;
    for (const auto& t : tracks)
    {
      if (!t.used)
      {
        continue;
      }
      const Acts::SquareMatrix2 nnT = t.wD0 * t.norm * t.norm.transpose();
      A += nnT;
      b += nnT * t.pos;
      ++n;
    }
    if (n < 2)
    {
      return false;
    }

    Eigen::LDLT<Acts::SquareMatrix2> ldlt(A);
    if (ldlt.info() != Eigen::Success || std::abs(A.determinant()) <= 0)
    {
      return false;
    }
    fit.xy = ldlt.solve(b);
    fit.covXY = A.inverse();

    double sumW = 0;
    double sumWZ = 0;
    for (const auto& t : tracks)
    {
      if (!t.used)
      {
        continue;
      }
      sumW += t.wZ0;
      sumWZ += t.wZ0 * t.zAt(fit.xy);
    }
    if (sumW <= 0)
    {
      return false;
    }
    fit.z = sumWZ / sumW;
    fit.varZ = 1. / sumW;
    fit.nTracks = n;
    return fit.xy.allFinite() && std::isfinite(fit.z);
  }
}  // namespace

//____________________________________________________________________________
ActsXYVertexSeeder::ActsXYVertexSeeder(Config cfg,
                                       std::unique_ptr<const Acts::Logger> logger)
  : m_cfg(std::move(cfg))
  , m_logger(std::move(logger))
{
  if (!m_cfg.zSeeder)
  {
    throw std::invalid_argument("ActsXYVertexSeeder: no z seed finder provided");
  }
  if (!m_cfg.extractParameters.connected())
  {
    throw std::invalid_argument("ActsXYVertexSeeder: no parameter extractor provided");
  }
}

//____________________________________________________________________________
Acts::Result<std::vector<Acts::Vertex>> ActsXYVertexSeeder::find(
    const std::vector<Acts::InputTrack>& trackVector,
    const Acts::VertexingOptions& vertexingOptions,
    Acts::IVertexFinder::State& /*state*/) const
{
  /// z estimate from the inner seeder. This is also the fallback seed
  auto zState = m_cfg.zSeeder->makeState(vertexingOptions.magFieldContext);
  auto zResult = m_cfg.zSeeder->find(trackVector, vertexingOptions, zState);
  if (!zResult.ok())
  {
    return zResult.error();
  }
  if (zResult->empty())
  {
    return std::vector<Acts::Vertex>{};
  }
  const Acts::Vertex fallback = zResult->back();

  /// beam spot as a weak prior in x,y, if its width is set
  const Acts::Vertex& constraint = vertexingOptions.constraint;
  const Acts::Vector2 priorXY = constraint.position().head<2>();
  Acts::SquareMatrix2 priorWeight = Acts::SquareMatrix2::Zero();
  const Acts::SquareMatrix2 priorCov = constraint.covariance().block<2, 2>(0, 0);
  if (priorCov(0, 0) > 0 && priorCov(1, 1) > 0)
  {
    priorWeight = priorCov.inverse();
  }

  /// track quantities
  std::vector<SeedTrack> tracks;
  tracks.reserve(trackVector.size());
  for (const auto& input : trackVector)
  {
    const Acts::BoundTrackParameters params = m_cfg.extractParameters(input);
    if (!params.covariance())
    {
      continue;
    }
    const auto& cov = *params.covariance();
    const double varD0 = cov(Acts::eBoundLoc0, Acts::eBoundLoc0);
    const double varZ0 = cov(Acts::eBoundLoc1, Acts::eBoundLoc1);
    const Acts::Vector3 direction = params.direction();
    const double dirT = direction.head<2>().norm();
    if (!(varD0 > 0) || !(varZ0 > 0) || dirT < 1e-6)
    {
      continue;
    }

    SeedTrack t;
    const Acts::Vector3 position = params.position(vertexingOptions.geoContext);
    t.pos = position.head<2>();
    t.z = position.z();
    t.dir = direction.head<2>() / dirT;
    t.norm = Acts::Vector2(-t.dir.y(), t.dir.x());
    t.cotTheta = direction.z() / dirT;
    t.wD0 = 1. / varD0;
    t.wZ0 = 1. / varZ0;
    tracks.push_back(t);
  }

  /// iterate: select tracks in the z window around the current seed, fit,
  /// drop the worst outlier and refit until all tracks are compatible
  SeedFit fit;
  fit.xy = priorXY;
  fit.z = fallback.position().z();
  bool haveFit = false;
  for (int iter = 0; iter < m_cfg.nIterations; ++iter)
  {
    for (auto& t : tracks)
    {
      t.used = std::abs(t.zAt(fit.xy) - fit.z) < m_cfg.zWindow;
    }

    SeedFit trial;
    bool ok = false;
    while (fitSeed(tracks, priorXY, priorWeight, trial))
    {
      double worstChi2 = -1;
      SeedTrack* worst = nullptr;
      for (auto& t : tracks)
      {
        if (!t.used)
        {
          continue;
        }
        const double chi2 = t.chi2(trial.xy, trial.z);
        if (chi2 > worstChi2)
        {
          worstChi2 = chi2;
          worst = &t;
        }
      }
      if (worstChi2 <= m_cfg.maxTrackChi2)
      {
        ok = true;
        break;
      }
      worst->used = false;
    }

    if (!ok)
    {
      // keep the result of the previous pass, if any
      break;
    }
    fit = trial;
    haveFit = true;

    ACTS_DEBUG("iteration " << iter << ": seed at (" << fit.xy.transpose()
                            << ", " << fit.z << ") from " << fit.nTracks
                            << " tracks");
  }

  if (!haveFit)
  {
    ACTS_DEBUG("transverse seed fit failed, using the z seed at "
               << fallback.fullPosition().transpose());
    return std::vector<Acts::Vertex>{fallback};
  }

  Acts::Vertex seed(Acts::Vector4(fit.xy.x(), fit.xy.y(), fit.z,
                                  fallback.fullPosition()[Acts::eTime]));
  Acts::SquareMatrix4 cov = Acts::SquareMatrix4::Zero();
  cov.block<2, 2>(0, 0) = fit.covXY;
  cov(2, 2) = fit.varZ;
  cov(3, 3) = constraint.fullCovariance()(3, 3);
  seed.setFullCovariance(cov);

  ACTS_DEBUG("seed at " << seed.fullPosition().transpose() << " from "
                        << fit.nTracks << " of " << trackVector.size()
                        << " tracks (z seed " << fallback.position().z() << ")");

  return std::vector<Acts::Vertex>{seed};
}
