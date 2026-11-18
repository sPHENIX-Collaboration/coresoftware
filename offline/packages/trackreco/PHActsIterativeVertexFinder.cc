#include "PHActsIterativeVertexFinder.h"

#include "ActsXYVertexSeeder.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrackVertexCrossingAssoc_v1.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <globalvertex/SvtxVertexMap_v1.h>
#include <globalvertex/SvtxVertex_v3.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/ParticleHypothesis.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Propagator.hpp>
#pragma GCC diagnostic pop
#include <Acts/Propagator/VoidNavigator.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Vertexing/FullBilloirVertexFitter.hpp>
#include <Acts/Vertexing/GaussianTrackDensity.hpp>
#include <Acts/Vertexing/HelicalTrackLinearizer.hpp>
#include <Acts/Vertexing/ImpactPointEstimator.hpp>
#include <Acts/Vertexing/IterativeVertexFinder.hpp>
#include <Acts/Vertexing/TrackAtVertex.hpp>
#include <Acts/Vertexing/TrackDensityVertexFinder.hpp>
#include <Acts/Vertexing/Vertex.hpp>
#include <Acts/Vertexing/VertexingOptions.hpp>
#include <Acts/Vertexing/ZScanVertexFinder.hpp>

#include <algorithm>
#include <climits>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>

namespace
{
  using Stepper = Acts::EigenStepper<>;
  using Propagator = Acts::Propagator<Stepper>;  // VoidNavigator by default

  /// time assigned to the tracks by ActsPropagator::makeTrackParams. Kept
  /// identical here so that the (unused) time dimension is consistent
  const double s_trackTime = 10 * Acts::UnitConstants::ns;

  /// tolerance for the on-surface check. The track position is placed on its
  /// own perigee line, so this only needs to absorb rounding
  const double s_onSurfaceTolerance = 1e-3 * Acts::UnitConstants::mm;

  Acts::Logging::Level actsLogLevel(int verbosity)
  {
    if (verbosity > 4)
    {
      return Acts::Logging::VERBOSE;
    }
    if (verbosity > 3)
    {
      return Acts::Logging::DEBUG;
    }
    if (verbosity > 0)
    {
      return Acts::Logging::WARNING;
    }
    return Acts::Logging::FATAL;
  }
}  // namespace

//____________________________________________________________________________
PHActsIterativeVertexFinder::PHActsIterativeVertexFinder(const std::string& name)
  : SubsysReco(name)
{
}

// out of line so that the unique_ptrs see the complete Acts types
PHActsIterativeVertexFinder::~PHActsIterativeVertexFinder() = default;

//____________________________________________________________________________
int PHActsIterativeVertexFinder::InitRun(PHCompositeNode* topNode)
{
  int ret = getNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  ret = createNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  setupFinder();

  if (Verbosity() > 0)
  {
    std::cout << PHWHERE << " configured with"
              << " seeder " << (m_seedFinder == SeedFinder::ZScan ? "ZScan" : m_seedFinder == SeedFinder::TrackDensity ? "TrackDensity" : "ZScanXYFit")
              << ", zScan minPt " << m_zScanMinPt << " GeV"
              << ", maxVertices " << m_maxVertices
              << ", significanceCutSeeding " << m_significanceCutSeeding
              << ", maxChi2CutForSeeding " << m_maxChi2CutForSeeding
              << ", reassign " << m_reassignTracksAfterFirstFit
              << ", beam spot (" << m_beamSpot.transpose() << ") cm"
              << ", beam spot width (" << m_beamSpotWidth.transpose() << ") cm"
              << ", beam spot constraint " << m_useBeamSpotConstraint
              << ", pp mode " << m_ppMode
              << ", outlier chi2 cut " << m_outlierChi2Cut
              << ", track cov scale " << m_trackCovScale
              << ", track MS term d0/z0 " << m_trackMSTermD0 << "/" << m_trackMSTermZ0 << " cm GeV"
              << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________
void PHActsIterativeVertexFinder::setupFinder()
{
  /// magnetic field: field map from the geometry, or a constant field
  m_field = m_tGeometry->geometry().magField;
  std::istringstream stringline(m_fieldMap);
  double fieldstrength = std::numeric_limits<double>::quiet_NaN();
  stringline >> fieldstrength;
  if (!m_fieldMap.empty() && !stringline.fail())
  {
    if (Verbosity() > 0)
    {
      std::cout << PHWHERE << " using constant field of " << fieldstrength << " T" << std::endl;
    }
    m_field = std::make_shared<Acts::ConstantBField>(
        Acts::Vector3(0, 0, fieldstrength * Acts::UnitConstants::T));
  }

  const auto logLevel = actsLogLevel(Verbosity());

  /// propagator without navigation, as in the Acts examples. Only used to
  /// transport track parameters between perigee surfaces near the beam line
  auto propagator = std::make_shared<Propagator>(
      Stepper(m_field), Acts::VoidNavigator{},
      Acts::getDefaultLogger("PHActsIVF::Propagator", logLevel));
  m_propagator = propagator;

  /// track linearizer
  Acts::HelicalTrackLinearizer::Config linearizerCfg;
  linearizerCfg.bField = m_field;
  linearizerCfg.propagator = m_propagator;
  m_linearizer = std::make_unique<Acts::HelicalTrackLinearizer>(
      linearizerCfg, Acts::getDefaultLogger("PHActsIVF::Linearizer", logLevel));

  /// vertex fitter
  Acts::FullBilloirVertexFitter::Config fitterCfg;
  fitterCfg.extractParameters.connect<&Acts::InputTrack::extractParameters>();
  fitterCfg.trackLinearizer.connect<&Acts::HelicalTrackLinearizer::linearizeTrack>(m_linearizer.get());
  Acts::FullBilloirVertexFitter fitter(
      fitterCfg, Acts::getDefaultLogger("PHActsIVF::Fitter", logLevel));
  /// identical fitter kept for the refits of the outlier pass
  m_fitter = std::make_unique<Acts::FullBilloirVertexFitter>(
      fitterCfg, Acts::getDefaultLogger("PHActsIVF::Refitter", logLevel));

  /// impact point estimator
  Acts::ImpactPointEstimator::Config ipEstCfg(m_field, m_propagator);
  Acts::ImpactPointEstimator ipEst(
      ipEstCfg, Acts::getDefaultLogger("PHActsIVF::IPEstimator", logLevel));

  /// seed finder
  std::shared_ptr<Acts::IVertexFinder> seeder;
  if (m_seedFinder == SeedFinder::ZScan || m_seedFinder == SeedFinder::ZScanXYFit)
  {
    Acts::ZScanVertexFinder::Config seederCfg(ipEst);
    seederCfg.extractParameters.connect<&Acts::InputTrack::extractParameters>();
    seederCfg.minPt = m_zScanMinPt * Acts::UnitConstants::GeV;
    seeder = std::make_shared<Acts::ZScanVertexFinder>(
        seederCfg, Acts::getDefaultLogger("PHActsIVF::ZScanSeeder", logLevel));

    if (m_seedFinder == SeedFinder::ZScanXYFit)
    {
      /// ZScan provides z, the x,y(,z) seed is fitted from the tracks
      ActsXYVertexSeeder::Config xyCfg;
      xyCfg.zSeeder = seeder;
      xyCfg.extractParameters.connect<&Acts::InputTrack::extractParameters>();
      xyCfg.zWindow = m_seedZWindow * Acts::UnitConstants::cm;
      xyCfg.maxTrackChi2 = m_seedMaxTrackChi2;
      seeder = std::make_shared<ActsXYVertexSeeder>(
          xyCfg, Acts::getDefaultLogger("PHActsIVF::XYSeeder", logLevel));
    }
  }
  else
  {
    Acts::GaussianTrackDensity::Config densityCfg;
    densityCfg.extractParameters.connect<&Acts::InputTrack::extractParameters>();
    seeder = std::make_shared<Acts::TrackDensityVertexFinder>(
        Acts::TrackDensityVertexFinder::Config{Acts::GaussianTrackDensity(densityCfg)});
  }

  /// the iterative vertex finder itself
  Acts::IterativeVertexFinder::Config finderCfg(std::move(fitter), seeder, ipEst);
  finderCfg.trackLinearizer.connect<&Acts::HelicalTrackLinearizer::linearizeTrack>(m_linearizer.get());
  finderCfg.extractParameters.connect<&Acts::InputTrack::extractParameters>();
  finderCfg.field = m_field;
  finderCfg.maxVertices = m_maxVertices;
  finderCfg.significanceCutSeeding = m_significanceCutSeeding;
  finderCfg.maximumChi2cutForSeeding = m_maxChi2CutForSeeding;
  finderCfg.reassignTracksAfterFirstFit = m_reassignTracksAfterFirstFit;
  finderCfg.cutOffTrackWeight = m_cutOffTrackWeight;

  m_finder = std::make_unique<Acts::IterativeVertexFinder>(
      std::move(finderCfg), Acts::getDefaultLogger("PHActsIVF", logLevel));
}

//____________________________________________________________________________
int PHActsIterativeVertexFinder::process_event(PHCompositeNode* /*topNode*/)
{
  ++m_nEvents;

  // in case these objects are in the input file, we clear the nodes and replace them
  m_vertexMap->Reset();
  m_trackVertexCrossingMap->Reset();

  /// sort tracks by crossing. Vertex ids are reset since the vertex map is
  std::set<short int> crossings;
  for (const auto& [trackkey, track] : *m_trackMap)
  {
    if (!track)
    {
      continue;
    }
    track->set_vertex_id(std::numeric_limits<unsigned int>::max());

    const short int crossing = track->get_crossing();

    // crossing zero contains unmatched TPC tracks
    // Here we skip those crossing = zero tracks that do not have silicon seeds
    if (m_ppMode && crossing == 0 && !getSiliconSeed(track))
    {
      continue;
    }

    crossings.insert(crossing);
    m_trackVertexCrossingMap->addTrackAssoc(crossing, trackkey);
  }

  unsigned int vertexId = 0;
  for (const auto crossing : crossings)
  {
    vertexId += processCrossing(crossing, vertexId);
  }

  /// update the crossing vertex map with the results
  for (const auto& [vtxkey, vertex] : *m_vertexMap)
  {
    m_trackVertexCrossingMap->addVertexAssoc(vertex->get_beam_crossing(), vtxkey);
  }

  if (Verbosity() > 0)
  {
    std::cout << PHWHERE << " event " << m_nEvents << ": " << m_trackMap->size()
              << " tracks, " << crossings.size() << " crossings, "
              << m_vertexMap->size() << " vertices" << std::endl;
  }
  if (Verbosity() > 1)
  {
    for (const auto& [vtxkey, vertex] : *m_vertexMap)
    {
      vertex->identify();
    }
  }
  if (Verbosity() > 2)
  {
    m_trackVertexCrossingMap->identify();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________
unsigned int PHActsIterativeVertexFinder::processCrossing(short int crossing,
                                                          unsigned int vertexId)
{
  const auto& geoContext = m_tGeometry->geometry().getGeoContext();
  const auto& magFieldContext = m_tGeometry->geometry().magFieldContext;

  /// collect the bound parameters at the fitted PCA for the selected
  /// tracks. The vector is filled completely before any pointer is taken so
  /// that the Acts::InputTrack pointers stay valid
  std::vector<Acts::BoundTrackParameters> trackParams;
  std::vector<unsigned int> trackKeys;

  ActsTransformations transformer;
  const auto range = m_trackVertexCrossingMap->getTracks(crossing);
  for (auto iter = range.first; iter != range.second; ++iter)
  {
    const unsigned int trackkey = iter->second;
    SvtxTrack* track = m_trackMap->get(trackkey);
    if (!track || !passTrackSelection(track))
    {
      continue;
    }

    SvtxTrackState* state = track->get_state(0.0);
    if (!state)
    {
      continue;
    }

    const Acts::Vector3 momentum(state->get_px(), state->get_py(), state->get_pz());
    const double p = momentum.norm();
    if (!std::isfinite(p) || p <= 0)
    {
      continue;
    }

    const Acts::Vector4 pos4(state->get_x() * Acts::UnitConstants::cm,
                             state->get_y() * Acts::UnitConstants::cm,
                             state->get_z() * Acts::UnitConstants::cm,
                             s_trackTime);

    /// round trip of PHActsTrkFitter's rotateActsCovToSvtxTrack, i.e. this
    /// recovers the fitted perigee covariance
    Acts::BoundSquareMatrix cov = transformer.rotateSvtxTrackCovToActs(state);

    /// optional diagnostic inflation of the track errors
    cov *= m_trackCovScale * m_trackCovScale;
    const double pt = std::hypot(state->get_px(), state->get_py());
    if (pt > 0)
    {
      const double bD0 = m_trackMSTermD0 * Acts::UnitConstants::cm / pt;
      const double bZ0 = m_trackMSTermZ0 * Acts::UnitConstants::cm / pt;
      cov(Acts::eBoundLoc0, Acts::eBoundLoc0) += bD0 * bD0;
      cov(Acts::eBoundLoc1, Acts::eBoundLoc1) += bZ0 * bZ0;
    }

    /// the vertex fitter inverts the covariance, reject ill-defined ones
    if (!cov.allFinite() || cov.llt().info() != Eigen::Success)
    {
      if (Verbosity() > 1)
      {
        std::cout << PHWHERE << " track " << trackkey
                  << " has a non positive definite covariance, skipping" << std::endl;
      }
      continue;
    }

    /// PHActsTrkFitter's reference surface is a perigee (a line along z)
    /// through the silicon seed's PCA to the origin, not through the origin
    /// itself, and that surface is not stored. A line along z through the
    /// track position is parallel to it and has the same PCA, so the bound
    /// covariance is identical and the position is on the surface by
    /// construction (d0 = z0 = 0). The IVF propagates from here to its own
    /// vertex candidates, so the choice of starting line does not matter
    const auto trackSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
        Acts::Vector3(pos4.head<3>()));

    auto result = Acts::BoundTrackParameters::create(
        geoContext, trackSurface, pos4, momentum / p,
        track->get_charge() / p, cov, Acts::ParticleHypothesis::pion(),
        s_onSurfaceTolerance);
    if (!result.ok())
    {
      if (Verbosity() > 1)
      {
        std::cout << PHWHERE << " could not create perigee parameters for track "
                  << trackkey << ": " << result.error().message() << std::endl;
      }
      continue;
    }

    if (Verbosity() > 3)
    {
      /// full input dump (Acts units: mm, GeV) so that a crossing can be
      /// replayed offline with exactly the same track parameters
      std::ostringstream dump;
      dump << std::setprecision(17) << "IVFINPUT crossing " << crossing << " key " << trackkey
           << " q " << track->get_charge()
           << " pos " << pos4.x() << " " << pos4.y() << " " << pos4.z()
           << " mom " << momentum.x() << " " << momentum.y() << " " << momentum.z()
           << " cov";
      for (int i = 0; i < 6; ++i)
      {
        for (int j = 0; j < 6; ++j)
        {
          dump << " " << cov(i, j);
        }
      }
      std::cout << dump.str() << std::endl;
    }

    trackParams.push_back(std::move(result.value()));
    trackKeys.push_back(trackkey);
  }

  if (Verbosity() > 0)
  {
    std::cout << PHWHERE << " crossing " << crossing << ": "
              << std::distance(range.first, range.second) << " tracks, "
              << trackParams.size() << " pass selection" << std::endl;
  }

  std::vector<unsigned int> vertexIds;

  /// the IVF needs at least two tracks to do anything
  if (trackParams.size() >= 2)
  {
    std::vector<Acts::InputTrack> inputTracks;
    std::map<const Acts::BoundTrackParameters*, unsigned int> paramsToKey;
    std::map<unsigned int, const Acts::BoundTrackParameters*> keyToParams;
    inputTracks.reserve(trackParams.size());
    for (std::size_t i = 0; i < trackParams.size(); ++i)
    {
      inputTracks.emplace_back(&trackParams[i]);
      paramsToKey.emplace(&trackParams[i], trackKeys[i]);
      keyToParams.emplace(trackKeys[i], &trackParams[i]);
    }

    /// beam spot. Its transverse position is the seed vertex position; it
    /// only enters the fit if the constraint is requested
    Acts::Vertex beamSpot(Acts::Vector4(m_beamSpot.x() * Acts::UnitConstants::cm,
                                        m_beamSpot.y() * Acts::UnitConstants::cm,
                                        m_beamSpot.z() * Acts::UnitConstants::cm,
                                        s_trackTime));
    Acts::SquareMatrix4 beamSpotCov = Acts::SquareMatrix4::Zero();
    beamSpotCov(0, 0) = std::pow(m_beamSpotWidth.x() * Acts::UnitConstants::cm, 2);
    beamSpotCov(1, 1) = std::pow(m_beamSpotWidth.y() * Acts::UnitConstants::cm, 2);
    beamSpotCov(2, 2) = std::pow(m_beamSpotWidth.z() * Acts::UnitConstants::cm, 2);
    // effectively unconstrained in time; must be non-zero as the fitter inverts it
    beamSpotCov(3, 3) = std::pow(1e6 * Acts::UnitConstants::ns, 2);
    beamSpot.setFullCovariance(beamSpotCov);

    const Acts::VertexingOptions options(geoContext, magFieldContext,
                                         beamSpot, m_useBeamSpotConstraint);

    auto state = m_finder->makeState(magFieldContext);

    ++m_nFinderCalls;
    auto result = m_finder->find(inputTracks, options, state);

    if (!result.ok())
    {
      ++m_nFinderFailures;
      if (Verbosity() > 0)
      {
        std::cout << PHWHERE << " Acts IVF failed for crossing " << crossing
                  << " with " << inputTracks.size() << " tracks: "
                  << result.error().message() << std::endl;
      }
    }
    else
    {
      if (result.value().empty() && Verbosity() > 1)
      {
        std::cout << PHWHERE << " crossing " << crossing
                  << ": Acts IVF returned no vertices from "
                  << inputTracks.size() << " tracks" << std::endl;
      }
      auto fieldCache = m_field->makeCache(magFieldContext);
      const bool doOutlierPass = m_outlierChi2Cut < std::numeric_limits<double>::max();

      for (const auto& vertex : result.value())
      {
        /// tracks attached to the vertex, each track used at most once
        std::vector<unsigned int> vtxTracks;
        unsigned int nLowWeight = 0;
        unsigned int nUnmatched = 0;
        unsigned int nAlreadyUsed = 0;
        for (const auto& trackAtVertex : vertex.tracks())
        {
          if (trackAtVertex.trackWeight < m_cutOffTrackWeight)
          {
            ++nLowWeight;
            continue;
          }
          const auto* params = trackAtVertex.originalParams.as<Acts::BoundTrackParameters>();
          const auto keyIter = paramsToKey.find(params);
          if (keyIter == paramsToKey.end())
          {
            ++nUnmatched;
            continue;
          }
          SvtxTrack* track = m_trackMap->get(keyIter->second);
          if (track->get_vertex_id() != std::numeric_limits<unsigned int>::max())
          {
            ++nAlreadyUsed;
            continue;
          }
          vtxTracks.push_back(keyIter->second);
        }

        if (vtxTracks.size() < m_minTracksPerVertex)
        {
          if (Verbosity() > 1)
          {
            std::cout << PHWHERE << " crossing " << crossing
                      << " dropping vertex at ("
                      << (vertex.position() / Acts::UnitConstants::cm).transpose()
                      << ") cm with chi2/ndf " << vertex.fitQuality().first
                      << "/" << vertex.fitQuality().second
                      << ": " << vtxTracks.size() << " usable of "
                      << vertex.tracks().size() << " attached tracks ("
                      << nLowWeight << " below weight " << m_cutOffTrackWeight
                      << ", " << nAlreadyUsed << " already used, "
                      << nUnmatched << " unmatched), need "
                      << m_minTracksPerVertex << std::endl;
          }
          continue;
        }

        /// outlier pass: remove the worst track and refit while its chi2 is
        /// above the cut. If tracks were dropped above (weight, already used)
        /// the vertex is first refitted with the tracks that are kept
        Acts::Vertex fitted = vertex;
        if (doOutlierPass)
        {
          bool needRefit = vtxTracks.size() != vertex.tracks().size();
          bool rejected = false;
          std::vector<unsigned int> outliers;
          while (true)
          {
            if (needRefit)
            {
              std::vector<Acts::InputTrack> fitTracks;
              fitTracks.reserve(vtxTracks.size());
              for (const auto key : vtxTracks)
              {
                fitTracks.emplace_back(keyToParams.at(key));
              }
              ++m_nRefits;
              auto refit = m_fitter->fit(fitTracks, options, fieldCache);
              if (!refit.ok())
              {
                ++m_nRefitFailures;
                if (Verbosity() > 1)
                {
                  std::cout << PHWHERE << " crossing " << crossing
                            << " refit with " << fitTracks.size() << " tracks failed: "
                            << refit.error().message() << std::endl;
                }
                rejected = true;
                break;
              }
              fitted = std::move(refit.value());
              needRefit = false;
            }

            /// largest track chi2 of the current fit decides whether a track
            /// has to go
            double worstChi2 = -1;
            for (const auto& trackAtVertex : fitted.tracks())
            {
              worstChi2 = std::max(worstChi2, trackAtVertex.chi2Track);
            }

            if (worstChi2 <= m_outlierChi2Cut)
            {
              break;
            }

            if (vtxTracks.size() <= m_minTracksPerVertex)
            {
              /// cannot remove another track
              rejected = true;
              break;
            }

            /// Which track goes is decided by leave-one-out: refit without
            /// each track in turn and remove the one whose removal gives the
            /// lowest vertex chi2. The track with the largest chi2 in the
            /// current fit is often a good one, because a single outlier
            /// with tiny errors drags the vertex away from all other tracks
            std::size_t bestIndex = vtxTracks.size();
            double bestChi2 = std::numeric_limits<double>::max();
            std::optional<Acts::Vertex> bestFit;
            for (std::size_t skip = 0; skip < vtxTracks.size(); ++skip)
            {
              std::vector<Acts::InputTrack> fitTracks;
              fitTracks.reserve(vtxTracks.size() - 1);
              for (std::size_t i = 0; i < vtxTracks.size(); ++i)
              {
                if (i != skip)
                {
                  fitTracks.emplace_back(keyToParams.at(vtxTracks[i]));
                }
              }
              ++m_nRefits;
              auto refit = m_fitter->fit(fitTracks, options, fieldCache);
              if (!refit.ok())
              {
                ++m_nRefitFailures;
                continue;
              }
              if (refit.value().fitQuality().first < bestChi2)
              {
                bestChi2 = refit.value().fitQuality().first;
                bestIndex = skip;
                bestFit = std::move(refit.value());
              }
            }

            if (!bestFit)
            {
              rejected = true;
              break;
            }

            const unsigned int removedKey = vtxTracks[bestIndex];
            if (Verbosity() > 1)
            {
              std::cout << PHWHERE << " crossing " << crossing << " vertex at ("
                        << (fitted.position() / Acts::UnitConstants::cm).transpose()
                        << ") cm, chi2 " << fitted.fitQuality().first
                        << ": removing outlier track " << removedKey
                        << ", chi2 without it " << bestChi2 << ", "
                        << vtxTracks.size() - 1 << " tracks left" << std::endl;
            }
            vtxTracks.erase(vtxTracks.begin() + static_cast<std::ptrdiff_t>(bestIndex));
            outliers.push_back(removedKey);
            ++m_nOutliersRemoved;
            fitted = std::move(*bestFit);
          }

          if (rejected)
          {
            ++m_nVerticesRejectedOutliers;
            if (Verbosity() > 1)
            {
              std::cout << PHWHERE << " crossing " << crossing
                        << " rejecting vertex at ("
                        << (fitted.position() / Acts::UnitConstants::cm).transpose()
                        << ") cm in outlier pass: " << vtxTracks.size()
                        << " tracks left after removing " << outliers.size()
                        << ", chi2/ndf " << fitted.fitQuality().first << "/"
                        << fitted.fitQuality().second << std::endl;
            }
            continue;
          }
        }

        const auto& [chi2, ndf] = fitted.fitQuality();
        if (ndf <= 0 || chi2 / ndf > m_maxVertexChi2Ndf)
        {
          if (Verbosity() > 1)
          {
            std::cout << PHWHERE << " rejecting vertex at ("
                      << (fitted.position() / Acts::UnitConstants::cm).transpose()
                      << ") cm with chi2/ndf " << chi2 << "/" << ndf << std::endl;
          }
          continue;
        }

        const unsigned int thisid = vertexId + vertexIds.size();

        auto svtxVertex = std::make_unique<SvtxVertex_v3>();
        svtxVertex->set_id(thisid);
        svtxVertex->set_beam_crossing(crossing);
        svtxVertex->set_t(0);

        const Acts::Vector3 position = fitted.position() / Acts::UnitConstants::cm;
        svtxVertex->set_x(position.x());
        svtxVertex->set_y(position.y());
        svtxVertex->set_z(position.z());

        const Acts::SquareMatrix3 vtxCov = fitted.covariance() / (Acts::UnitConstants::cm * Acts::UnitConstants::cm);
        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
          {
            svtxVertex->set_error(i, j, vtxCov(i, j));
          }
        }

        /// note that PHActsVertexPropagator::setVtxChi2 overwrites these
        svtxVertex->set_chisq(chi2);
        svtxVertex->set_ndof(ndf > 0 ? static_cast<unsigned int>(std::lround(ndf)) : 0);

        for (const auto trackkey : vtxTracks)
        {
          svtxVertex->insert_track(trackkey);
          m_trackMap->get(trackkey)->set_vertex_id(thisid);
        }

        if (Verbosity() > 1)
        {
          std::cout << PHWHERE << " crossing " << crossing << " vertex " << thisid
                    << " at (" << position.transpose() << ") cm with "
                    << vtxTracks.size() << " tracks, chi2/ndf " << chi2 << "/" << ndf
                    << std::endl;
        }
        if (Verbosity() > 2)
        {
          /// per-track Billoir chi2 contribution, to spot single outlier tracks
          for (const auto& trackAtVertex : fitted.tracks())
          {
            const auto* params = trackAtVertex.originalParams.as<Acts::BoundTrackParameters>();
            const auto keyIter = paramsToKey.find(params);
            std::cout << "    track " << (keyIter == paramsToKey.end() ? -1 : static_cast<long>(keyIter->second))
                      << " pT " << params->transverseMomentum()
                      << " sigma(d0,z0) (" << std::sqrt((*params->covariance())(Acts::eBoundLoc0, Acts::eBoundLoc0)) / Acts::UnitConstants::um
                      << ", " << std::sqrt((*params->covariance())(Acts::eBoundLoc1, Acts::eBoundLoc1)) / Acts::UnitConstants::um
                      << ") um chi2 " << trackAtVertex.chi2Track
                      << " weight " << trackAtVertex.trackWeight << std::endl;
          }
        }

        m_vertexMap->insert(svtxVertex.release());
        vertexIds.push_back(thisid);
      }
    }
  }

  m_nVertices += vertexIds.size();

  if (m_assignUnusedTracks && !vertexIds.empty())
  {
    assignUnusedTracks(crossing, vertexIds);
  }

  return vertexIds.size();
}

//____________________________________________________________________________
void PHActsIterativeVertexFinder::assignUnusedTracks(short int crossing,
                                                     const std::vector<unsigned int>& vertexIds)
{
  /// Assign the closest vertex in z to every track of this crossing that was
  /// not used by any vertex, so that PHActsVertexPropagator can propagate it.
  /// The track is not added to the vertex track list
  const auto range = m_trackVertexCrossingMap->getTracks(crossing);
  for (auto iter = range.first; iter != range.second; ++iter)
  {
    SvtxTrack* track = m_trackMap->get(iter->second);
    if (!track || track->get_vertex_id() != std::numeric_limits<unsigned int>::max())
    {
      continue;
    }

    float mindz = std::numeric_limits<float>::max();
    unsigned int newvtxid = std::numeric_limits<unsigned int>::max();
    for (const auto vtxid : vertexIds)
    {
      const float dz = std::abs(track->get_z() - m_vertexMap->get(vtxid)->get_z());
      if (dz < mindz)
      {
        mindz = dz;
        newvtxid = vtxid;
      }
    }

    track->set_vertex_id(newvtxid);
    if (Verbosity() > 2)
    {
      std::cout << PHWHERE << " assign vertex " << newvtxid
                << " to unused track " << iter->second << std::endl;
    }
  }
}

//____________________________________________________________________________
bool PHActsIterativeVertexFinder::passTrackSelection(SvtxTrack* track) const
{
  if (track->get_quality() > m_qualCut)
  {
    return false;
  }
  if (track->get_pt() < m_trackPtCut)
  {
    return false;
  }
  if (m_requireMvtx && countSiliconClusters(track, TrkrDefs::mvtxId) < m_nMvtxRequired)
  {
    return false;
  }
  if (m_requireIntt && countSiliconClusters(track, TrkrDefs::inttId) < m_nInttRequired)
  {
    return false;
  }
  return true;
}

//____________________________________________________________________________
unsigned int PHActsIterativeVertexFinder::countSiliconClusters(SvtxTrack* track, uint8_t trkrId) const
{
  TrackSeed* siliconseed = getSiliconSeed(track);
  if (!siliconseed)
  {
    return 0;
  }

  unsigned int nclus = 0;
  for (auto clusit = siliconseed->begin_cluster_keys(); clusit != siliconseed->end_cluster_keys(); ++clusit)
  {
    if (TrkrDefs::getTrkrId(*clusit) == trkrId)
    {
      ++nclus;
    }
  }
  return nclus;
}

//____________________________________________________________________________
TrackSeed* PHActsIterativeVertexFinder::getSiliconSeed(SvtxTrack* track) const
{
  // SvtxTrack_v4 and earlier carry the seed pointer, SvtxTrack_v5 only the seed id
  if (TrackSeed* seed = track->get_silicon_seed())
  {
    return seed;
  }
  const unsigned short id = track->get_silicon_seed_id();
  // find() does no bounds checking, so guard against unset (USHRT_MAX) or stale ids
  if (!m_siliconSeedMap || id >= m_siliconSeedMap->size())
  {
    return nullptr;
  }
  return *m_siliconSeedMap->find(id);
}

//____________________________________________________________________________
int PHActsIterativeVertexFinder::End(PHCompositeNode* /*topNode*/)
{
  std::cout << PHWHERE << " processed " << m_nEvents << " events, "
            << m_nFinderCalls << " IVF calls (" << m_nFinderFailures << " failed), "
            << m_nVertices << " vertices found" << std::endl;
  std::cout << PHWHERE << " outlier pass: " << m_nOutliersRemoved << " tracks removed, "
            << m_nRefits << " refits (" << m_nRefitFailures << " failed), "
            << m_nVerticesRejectedOutliers << " vertices rejected" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________
int PHActsIterativeVertexFinder::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << " DST node missing, quitting" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  auto* svtxNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, m_vertexMapName);
  if (!m_vertexMap)
  {
    m_vertexMap = new SvtxVertexMap_v1;
    auto* vertexNode = new PHIODataNode<PHObject>(m_vertexMap, m_vertexMapName, "PHObject");
    svtxNode->addNode(vertexNode);
  }

  m_trackVertexCrossingMap = findNode::getClass<TrackVertexCrossingAssoc>(topNode, "TrackVertexCrossingAssocMap");
  if (!m_trackVertexCrossingMap)
  {
    m_trackVertexCrossingMap = new TrackVertexCrossingAssoc_v1;
    auto* assocNode = new PHIODataNode<PHObject>(m_trackVertexCrossingMap, "TrackVertexCrossingAssocMap", "PHObject");
    svtxNode->addNode(assocNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________
int PHActsIterativeVertexFinder::getNodes(PHCompositeNode* topNode)
{
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  if (!m_trackMap)
  {
    std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << m_trackMapName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE << " ERROR: Can't find acts tracking geometry" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_siliconSeedMap = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if (!m_siliconSeedMap)
  {
    std::cout << PHWHERE << " ERROR: Can't find SiliconTrackSeedContainer" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
