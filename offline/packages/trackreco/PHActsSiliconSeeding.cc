#include "PHActsSiliconSeeding.h"

#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase_historic/ActsTransformations.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <intt/CylinderGeomIntt.h>

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterCrossingAssoc.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeedHelper.h>
#include <trackbase_historic/TrackSeed_v2.h>

#ifndef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif
#include <Acts/Seeding/BinnedSPGroup.hpp>
#ifndef __clang__
#pragma GCC diagnostic pop
#endif
#include <Acts/Seeding/InternalSeed.hpp>
#include <Acts/Seeding/InternalSpacePoint.hpp>
#include <Acts/Seeding/Seed.hpp>
#include <Acts/Seeding/SeedFilter.hpp>

namespace
{
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }
}  // namespace

PHActsSiliconSeeding::PHActsSiliconSeeding(const std::string& name)
  : SubsysReco(name)
{
}

int PHActsSiliconSeeding::Init(PHCompositeNode* /*topNode*/)
{
  Acts::SeedFilterConfig sfCfg = configureSeedFilter();
  sfCfg = sfCfg.toInternalUnits();

  m_seedFinderCfg.seedFilter = std::make_unique<Acts::SeedFilter<SpacePoint>>(
      Acts::SeedFilter<SpacePoint>(sfCfg));

  configureSeeder();
  configureSPGrid();

  if (Verbosity() > 5)
  {
    printSeedConfigs(sfCfg);
  }
  // vector containing the map of z bins in the top and bottom layers
  m_bottomBinFinder = std::make_shared<const Acts::BinFinder<SpacePoint>>(
      zBinNeighborsBottom, nphineighbors);
  m_topBinFinder = std::make_shared<const Acts::BinFinder<SpacePoint>>(
      zBinNeighborsTop, nphineighbors);

  if (m_seedAnalysis)
  {
    m_file = new TFile(std::string(Name() + ".root").c_str(), "recreate");
    createHistograms();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
int PHActsSiliconSeeding::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsSiliconSeeding::process_event(PHCompositeNode* topNode)
{
  if (m_nIteration > 0)
  {
    _iteration_map = findNode::getClass<TrkrClusterIterationMapv1>(topNode, "TrkrClusterIterationMap");
    if (!_iteration_map)
    {
      std::cerr << PHWHERE << "Cluster Iteration Map missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  auto eventTimer = std::make_unique<PHTimer>("eventTimer");
  eventTimer->stop();
  eventTimer->restart();

  if (Verbosity() > 0)
  {
    std::cout << "Processing PHActsSiliconSeeding event "
              << m_event << std::endl;
  }

  std::vector<const SpacePoint*> spVec;
  auto seedVector = runSeeder(spVec);

  eventTimer->stop();
  auto seederTime = eventTimer->get_accumulated_time();
  eventTimer->restart();

  makeSvtxTracks(seedVector);

  eventTimer->stop();
  auto circleFitTime = eventTimer->get_accumulated_time();

  for (auto sp : spVec)
  {
    delete sp;
  }
  spVec.clear();

  if (Verbosity() > 0)
  {
    std::cout << "Finished PHActsSiliconSeeding process_event"
              << std::endl;
  }

  if (Verbosity() > 0)
  {
    std::cout << "PHActsSiliconSeeding Acts seed time "
              << seederTime << std::endl;
    std::cout << "PHActsSiliconSeeding circle fit time "
              << circleFitTime << std::endl;
    std::cout << "PHActsSiliconSeeding total event time "
              << circleFitTime + seederTime << std::endl;
  }

  m_event++;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsSiliconSeeding::End(PHCompositeNode* /*topNode*/)
{
  if (m_seedAnalysis)
  {
    writeHistograms();
  }

  if (Verbosity() > 1)
  {
    std::cout << "There were " << m_nBadInitialFits
              << " bad initial circle fits" << std::endl;
    std::cout << "There were " << m_nBadUpdates
              << " bad second circle fits" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

GridSeeds PHActsSiliconSeeding::runSeeder(std::vector<const SpacePoint*>& spVec)
{
  Acts::SeedFinder<SpacePoint> seedFinder(m_seedFinderCfg);
  GridSeeds seedVector;
  for (int strobe = m_lowStrobeIndex; strobe < m_highStrobeIndex; strobe++)
  {
    /// Covariance converter functor needed by seed finder
    auto covConverter =
        [=](const SpacePoint& sp, float zAlign, float rAlign, float sigmaError)
        -> std::pair<Acts::Vector3, Acts::Vector2>
    {
      Acts::Vector3 position{sp.x(), sp.y(), sp.z()};
      Acts::Vector2 cov;
      cov[0] = (sp.m_varianceR + rAlign * rAlign) * sigmaError;
      cov[1] = (sp.m_varianceZ + zAlign * zAlign) * sigmaError;
      return std::make_pair(position, cov);
    };

    Acts::Extent rRangeSPExtent;

    spVec = getSiliconSpacePoints(rRangeSPExtent, strobe);

    if (m_seedAnalysis)
    {
      h_nInputMeas->Fill(spVec.size());
    }

    auto grid =
        Acts::SpacePointGridCreator::createGrid<SpacePoint>(m_gridCfg,
                                                            m_gridOptions);

    auto spGroup = Acts::BinnedSPGroup<SpacePoint>(spVec.begin(),
                                                   spVec.end(),
                                                   covConverter,
                                                   m_bottomBinFinder,
                                                   m_topBinFinder,
                                                   std::move(grid),
                                                   rRangeSPExtent,
                                                   m_seedFinderCfg,
                                                   m_seedFinderOptions);

    /// variable middle SP radial region of interest
    const Acts::Range1D<float> rMiddleSPRange(
        std::floor(rRangeSPExtent.min(Acts::binR) / 2) * 2 + 1.5,
        std::floor(rRangeSPExtent.max(Acts::binR) / 2) * 2 - 1.5);

    SeedContainer seeds;
    seeds.clear();
    decltype(seedFinder)::SeedingState state;
    state.spacePointData.resize(spVec.size(),
                                m_seedFinderCfg.useDetailedDoubleMeasurementInfo);
    for (const auto [bottom, middle, top] : spGroup)
    {
      seedFinder.createSeedsForGroup(m_seedFinderOptions,
                                     state, spGroup.grid(),
                                     std::back_inserter(seeds),
                                     bottom,
                                     middle,
                                     top,
                                     rMiddleSPRange);
    }

    seedVector.push_back(seeds);
  }

  return seedVector;
}

void PHActsSiliconSeeding::makeSvtxTracks(GridSeeds& seedVector)
{
  int numSeeds = 0;
  int numGoodSeeds = 0;
  m_seedid = -1;
  int strobe = m_lowStrobeIndex;
  /// Loop over grid volumes. In our case this will be strobe
  for (auto& seeds : seedVector)
  {
    /// Loop over actual seeds in this grid volume
    for (auto& seed : seeds)
    {
      if (Verbosity() > 1)
      {
        std::cout << "Seed " << numSeeds << " has "
                  << seed.sp().size() << " measurements "
                  << std::endl;
      }

      if (m_seedAnalysis)
      {
        clearTreeVariables();
        m_seedid++;
      }

      std::vector<TrkrDefs::cluskey> cluster_keys;

      std::vector<Acts::Vector3> globalPositions;

      std::map<TrkrDefs::cluskey, Acts::Vector3> positions;
      auto trackSeed = std::make_unique<TrackSeed_v2>();

      for (auto& spacePoint : seed.sp())
      {
        const auto& cluskey = spacePoint->Id();
        cluster_keys.push_back(cluskey);

        trackSeed->insert_cluster_key(cluskey);
        auto globalPosition = m_tGeometry->getGlobalPosition(
            cluskey,
            m_clusterMap->findCluster(cluskey));
        globalPositions.push_back(globalPosition);
        if (m_seedAnalysis)
        {
          m_mvtxgx.push_back(globalPosition(0));
          m_mvtxgy.push_back(globalPosition(1));
          float clusr = std::sqrt(square(globalPosition(0)) + square(globalPosition(1)));
          if (globalPosition.y() < 0)
          {
            clusr *= -1;
          }
          m_mvtxgr.push_back(clusr);
          m_mvtxgz.push_back(globalPosition(2));
        }
        positions.insert(std::make_pair(cluskey, globalPosition));
        if (Verbosity() > 1)
        {
          std::cout << "Adding cluster with x,y "
                    << spacePoint->x() << ", " << spacePoint->y()
                    << " mm in detector "
                    << (unsigned int) TrkrDefs::getTrkrId(cluskey)
                    << " with cluskey " << cluskey
                    << std::endl;
        }
      }
      if (m_searchInIntt)
      {
        int nintt = 0;
        for (auto& key : cluster_keys)
        {
          if (TrkrDefs::getTrkrId(key) == TrkrDefs::TrkrId::inttId)
          {
            nintt++;
          }
        }
        /// if acts found a triplet in the INTT only it is likely a bad seed
        if (nintt > 2)
        {
          continue;
        }
      }
      double z = seed.z() / Acts::UnitConstants::cm;

      auto fitTimer = std::make_unique<PHTimer>("trackfitTimer");
      fitTimer->stop();
      fitTimer->restart();

      TrackSeedHelper::circleFitByTaubin(trackSeed.get(), positions, 0, 8);
      const auto position(TrackSeedHelper::get_xyz(trackSeed.get()));
      if (std::abs(position.x()) > m_maxSeedPCA || std::abs(position.y()) > m_maxSeedPCA)
      {
        if (Verbosity() > 1)
        {
          std::cout << "Large PCA seed " << std::endl;
          trackSeed->identify();
        }
        m_nBadInitialFits++;
        continue;
      }

      TrackSeedHelper::lineFit(trackSeed.get(), positions, 0, 8);
      z = trackSeed->get_Z0();
      fitTimer->stop();
      auto circlefittime = fitTimer->get_accumulated_time();
      fitTimer->restart();

      // calculate phi and assign
      auto phi = TrackSeedHelper::get_phi(trackSeed.get(), positions);
      trackSeed->set_phi(phi);

      /// Project to INTT and find matches
      const auto mvtxsize = globalPositions.size();
      auto additionalClusters = findMatches(globalPositions, cluster_keys, *trackSeed);

      /// Add possible matches to cluster list to be parsed when
      /// Svtx tracks are made
      for (unsigned int newkey = 0; newkey < additionalClusters.size(); newkey++)
      {
        trackSeed->insert_cluster_key(additionalClusters[newkey]);
        positions.insert(std::make_pair(additionalClusters[newkey], globalPositions[mvtxsize + newkey]));

        if (Verbosity() > 2)
        {
          std::cout << "adding additional intt key " << additionalClusters[newkey] << std::endl;
        }
      }

      fitTimer->stop();
      auto addClusters = fitTimer->get_accumulated_time();
      fitTimer->restart();

      //! Circle fit again to take advantage of INTT lever arm
      TrackSeedHelper::circleFitByTaubin(trackSeed.get(), positions, 0, 7);
      phi = TrackSeedHelper::get_phi(trackSeed.get(), positions);
      trackSeed->set_phi(phi);
      if (m_searchInIntt)
      {
        TrackSeedHelper::lineFit(trackSeed.get(), positions, 0, 2);
      }

      if (Verbosity() > 0)
      {
        std::cout << "find intt clusters time " << addClusters << std::endl;
      }

      numGoodSeeds++;

      /// The Acts z projection has better resolution than the circle fit
      trackSeed->set_Z0(z);

      if (Verbosity() > 1)
      {
        std::cout << "Silicon seed id " << m_seedContainer->size() << std::endl;
        std::cout << "seed phi, theta, eta : "
                  << trackSeed->get_phi() << ", " << trackSeed->get_theta()
                  << ", " << trackSeed->get_eta() << std::endl;
        trackSeed->identify();
      }

      //! try to get a crossing value based on INTT
      trackSeed->set_crossing(getCrossingIntt(*trackSeed));

      m_seedContainer->insert(trackSeed.get());

      fitTimer->stop();
      auto svtxtracktime = fitTimer->get_accumulated_time();
      if (Verbosity() > 0)
      {
        std::cout << "Intt fit time " << circlefittime << " and svtx time "
                  << svtxtracktime << std::endl;
      }
    }
    strobe++;
    if (strobe > m_highStrobeIndex)
    {
      std::cout << PHWHERE << "Error: some how grid seed vector is not the same as the number of strobes" << std::endl;
    }
  }

  if (m_seedAnalysis)
  {
    h_nSeeds->Fill(numGoodSeeds);
    h_nActsSeeds->Fill(numSeeds);
  }
  if (Verbosity() > 1)
  {
    std::cout << "Total number of seeds found in "
              << seedVector.size() << " volume regions gives "
              << numSeeds << " Acts seeds " << std::endl
              << std::endl
              << std::endl;
    m_seedContainer->identify();
    for (auto& seed : *m_seedContainer)
    {
      if (!seed)
      {
        continue;
      }
      seed->identify();
    }
  }

  return;
}

short int PHActsSiliconSeeding::getCrossingIntt(TrackSeed& si_track)
{
  // If the Si track contains an INTT hit, use it to get the bunch crossing offset

  std::vector<short int> intt_crossings = getInttCrossings(si_track);

  bool keep_it = true;
  short int crossing_keep = 0;
  if (intt_crossings.size() == 0)
  {
    keep_it = false;
  }
  else
  {
    crossing_keep = intt_crossings[0];
    for (unsigned int ic = 1; ic < intt_crossings.size(); ++ic)
    {
      if (intt_crossings[ic] != crossing_keep)
      {
        if (Verbosity() > 1)
        {
          std::cout << " Warning: INTT crossings not all the same "
                    << " crossing_keep " << crossing_keep << " new crossing " << intt_crossings[ic] << " keep the first one in the list" << std::endl;
        }
      }
    }
  }

  if (keep_it)
  {
    return crossing_keep;
  }

  return SHRT_MAX;
}

std::vector<short int> PHActsSiliconSeeding::getInttCrossings(TrackSeed& si_track)
{
  std::vector<short int> intt_crossings;

  // If the Si track contains an INTT hit, use it to get the bunch crossing offset
  // loop over associated clusters to get keys for silicon cluster
  for (TrackSeed::ConstClusterKeyIter iter = si_track.begin_cluster_keys();
       iter != si_track.end_cluster_keys();
       ++iter)
  {
    TrkrDefs::cluskey cluster_key = *iter;
    const unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);

    if (Verbosity() > 1)
    {
      unsigned int layer = TrkrDefs::getLayer(cluster_key);

      if (trkrid == TrkrDefs::mvtxId)
      {
        TrkrCluster* cluster = m_clusterMap->findCluster(cluster_key);
        if (!cluster)
        {
          continue;
        }

        Acts::Vector3 global = m_tGeometry->getGlobalPosition(cluster_key, cluster);

        std::cout << "Checking  si Track with cluster " << cluster_key
                  << " in layer " << layer << " position " << global(0) << "  " << global(1) << "  " << global(2)
                  << " eta " << si_track.get_eta() << std::endl;
      }
      else
      {
        std::cout << "Checking  si Track with cluster " << cluster_key
                  << " in layer " << layer << " with eta " << si_track.get_eta() << std::endl;
      }
    }

    if (trkrid == TrkrDefs::inttId)
    {
      TrkrCluster* cluster = m_clusterMap->findCluster(cluster_key);
      if (!cluster)
      {
        continue;
      }

      unsigned int layer = TrkrDefs::getLayer(cluster_key);

      // get the bunch crossings for all hits in this cluster
      auto crossings = _cluster_crossing_map->getCrossings(cluster_key);
      for (auto iter1 = crossings.first; iter1 != crossings.second; ++iter1)
      {
        if (Verbosity() > 1)
        {
          std::cout << "                si Track with cluster " << iter1->first << " layer " << layer << " crossing " << iter1->second << std::endl;
        }
        intt_crossings.push_back(iter1->second);
      }
    }
  }

  return intt_crossings;
}

std::vector<TrkrDefs::cluskey> PHActsSiliconSeeding::findMatches(
    std::vector<Acts::Vector3>& clusters,
    std::vector<TrkrDefs::cluskey>& keys,
    TrackSeed& seed)
{
  auto fitpars = TrackFitUtils::fitClusters(clusters, keys, true);

  float trackphi = seed.get_phi();
  /// Diagnostic
  if (m_seedAnalysis)
  {
    for (const auto& glob : clusters)
    {
      h_hits->Fill(glob(0), glob(1));
      h_zhits->Fill(glob(2),
                    std::sqrt(square(glob(0)) + square(glob(1))));
      h_projHits->Fill(glob(0), glob(1));
      h_zprojHits->Fill(glob(2),
                        sqrt(square(glob(0)) + square(glob(1))));
    }
  }

  std::vector<TrkrDefs::cluskey> matchedClusters;
  std::map<int, float> minResidLayer;
  std::map<int, TrkrDefs::cluskey> minResidckey;
  std::map<int, Acts::Vector3> minResidGlobPos;
  for (int i = 0; i < 7; i++)
  {
    minResidLayer.insert(std::make_pair(i, std::numeric_limits<float>::max()));
  }
  std::set<unsigned int> layersToSkip;
  for (auto it = seed.begin_cluster_keys();
       it != seed.end_cluster_keys();
       ++it)
  {
    auto key = *it;
    unsigned int layer = TrkrDefs::getLayer(key);
    layersToSkip.insert(layer);
  }

  int nlayers = 3;
  int layer = 0;
  for (auto& det : {TrkrDefs::TrkrId::mvtxId, TrkrDefs::TrkrId::inttId})
  {
    if (det == TrkrDefs::TrkrId::inttId)
    {
      nlayers = 7;
      layer = 3;
    }
    while (layer < nlayers)
    {
      if (layersToSkip.find(layer) != layersToSkip.end())
      {
        layer++;
        continue;
      }
      for (const auto& hitsetkey : m_clusterMap->getHitSetKeys(det, layer))
      {
        auto surf = m_tGeometry->maps().getSiliconSurface(hitsetkey);
        auto surfcenter = surf->center(m_tGeometry->geometry().geoContext);
        float surfphi = atan2(surfcenter.y(), surfcenter.x());
        float dphi = normPhi2Pi(trackphi - surfphi);

        /// Check that the projection is within some reasonable amount of the segment
        /// to reject e.g. looking at segments in the opposite hemisphere. This is about
        /// the size of one intt segment (256 * 80 micron strips in a segment)
        if (fabs(dphi) > 0.2)
        {
          continue;
        }

        auto range = m_clusterMap->getClusters(hitsetkey);
        for (auto clusIter = range.first; clusIter != range.second; ++clusIter)
        {
          const auto cluskey = clusIter->first;
          if (_iteration_map != nullptr && m_nIteration > 0)
          {
            if (_iteration_map->getIteration(cluskey) < m_nIteration)
            {
              continue;  // skip clusters used in a previous iteration
            }
          }

          const auto cluster = clusIter->second;
          auto glob = m_tGeometry->getGlobalPosition(
              cluskey, cluster);
          auto intersection = TrackFitUtils::get_helix_surface_intersection(surf, fitpars, glob, m_tGeometry);
          auto local = (surf->transform(m_tGeometry->geometry().getGeoContext())).inverse() * (intersection * Acts::UnitConstants::cm);
          local /= Acts::UnitConstants::cm;
          m_projgx = intersection.x();
          m_projgy = intersection.y();
          m_projgr = std::sqrt(square(intersection.x()) + square(intersection.y()));
          if (m_projgy < 0)
          {
            m_projgr *= -1;
          }
          m_projgz = intersection.z();
          m_projlx = local.x();
          m_projlz = local.y();
          /// Diagnostic
          if (m_seedAnalysis)
          {
            const auto globalP = m_tGeometry->getGlobalPosition(
                cluskey, cluster);
            m_clusgx = globalP.x();
            m_clusgy = globalP.y();
            m_clusgr = std::sqrt(square(globalP.x()) + square(globalP.y()));
            if (globalP.y() < 0)
            {
              m_clusgr *= -1;
            }
            m_clusgz = globalP.z();
            m_cluslx = cluster->getLocalX();
            m_cluslz = cluster->getLocalY();
            h_nInttProj->Fill(local.x() - cluster->getLocalX(),
                              local.y() - cluster->getLocalY());
            h_hits->Fill(globalP(0), globalP(1));
            h_zhits->Fill(globalP(2),
                          std::sqrt(square(globalP(0)) + square(globalP(1))));

            h_resids->Fill(local.y() - cluster->getLocalY(),
                           local.x() - cluster->getLocalX());
            m_tree->Fill();
          }

          /// Z strip spacing is the entire strip, so because we use fabs
          /// we divide by two
          float rphiresid = fabs(local.x() - cluster->getLocalX());
          float zresid = fabs(local.y() - cluster->getLocalY());

          if ((det == TrkrDefs::TrkrId::mvtxId && rphiresid < m_mvtxrPhiSearchWin &&
               zresid < m_mvtxzSearchWin) ||
              (det == TrkrDefs::TrkrId::inttId && rphiresid < m_inttrPhiSearchWin && zresid < m_inttzSearchWin))

          {
            if (rphiresid < minResidLayer[layer])
            {
              minResidLayer[layer] = rphiresid;
              minResidckey[layer] = cluskey;
              minResidGlobPos[layer] = glob;
            }

            if (Verbosity() > 4)
            {
              std::cout << "Matched INTT cluster with cluskey " << cluskey
                        << " and position " << glob.transpose()
                        << std::endl
                        << " with projections rphi "
                        << local.x() << " and inttclus rphi " << cluster->getLocalX()
                        << " and proj z " << local.y() << " and inttclus z "
                        << cluster->getLocalY() << " in layer " << layer
                        << std::endl;
            }
          }
        }
      }
      layer++;
    }
  }
  for (int ilayer = 0; ilayer < 3; ilayer++)
  {
    if (minResidLayer[ilayer] < std::numeric_limits<float>::max())
    {
      matchedClusters.push_back(minResidckey[ilayer]);
      clusters.push_back(minResidGlobPos[ilayer]);
    }
  }

  if (minResidLayer[3] < minResidLayer[4] && minResidLayer[3] < std::numeric_limits<float>::max())
  {
    matchedClusters.push_back(minResidckey[3]);
    clusters.push_back(minResidGlobPos[3]);
  }
  else if (minResidLayer[4] < std::numeric_limits<float>::max())
  {
    matchedClusters.push_back(minResidckey[4]);
    clusters.push_back(minResidGlobPos[4]);
  }

  if (minResidLayer[5] < minResidLayer[6] && minResidLayer[5] < std::numeric_limits<float>::max())
  {
    matchedClusters.push_back(minResidckey[5]);
    clusters.push_back(minResidGlobPos[5]);
  }
  else if (minResidLayer[6] < std::numeric_limits<float>::max())
  {
    matchedClusters.push_back(minResidckey[6]);
    clusters.push_back(minResidGlobPos[6]);
  }

  if (m_seedAnalysis)
  {
    h_nMatchedClusters->Fill(matchedClusters.size());
  }

  return matchedClusters;
}

SpacePointPtr PHActsSiliconSeeding::makeSpacePoint(
    const Surface& surf,
    const TrkrDefs::cluskey key,
    TrkrCluster* clus)
{
  Acts::Vector2 localPos(clus->getLocalX() * Acts::UnitConstants::cm,
                         clus->getLocalY() * Acts::UnitConstants::cm);
  Acts::Vector3 globalPos(0, 0, 0);
  Acts::Vector3 mom(1, 1, 1);
  globalPos = surf->localToGlobal(m_tGeometry->geometry().getGeoContext(),
                                  localPos, mom);

  Acts::SquareMatrix2 localCov = Acts::SquareMatrix2::Zero();
  localCov(0, 0) = pow(clus->getRPhiError(), 2) * Acts::UnitConstants::cm2;
  localCov(1, 1) = pow(clus->getZError(), 2) * Acts::UnitConstants::cm2;

  float x = globalPos.x();
  float y = globalPos.y();
  float z = globalPos.z();
  float r = std::sqrt(x * x + y * y);

  /// The space point requires only the variance of the transverse and
  /// longitudinal position. Reduce computations by transforming the
  /// covariance directly from local to r/z.
  ///
  /// compute Jacobian from global coordinates to r/z
  ///
  ///         r = sqrt(x² + y²)
  /// dr/d{x,y} = (1 / sqrt(x² + y²)) * 2 * {x,y}
  ///             = 2 * {x,y} / r
  ///       dz/dz = 1

  Acts::RotationMatrix3 rotLocalToGlobal =
      surf->referenceFrame(m_tGeometry->geometry().getGeoContext(), globalPos, mom);
  auto scale = 2 / std::hypot(x, y);
  Acts::ActsMatrix<2, 3> jacXyzToRhoZ = Acts::ActsMatrix<2, 3>::Zero();
  jacXyzToRhoZ(0, Acts::ePos0) = scale * x;
  jacXyzToRhoZ(0, Acts::ePos1) = scale * y;
  jacXyzToRhoZ(1, Acts::ePos2) = 1;
  // compute Jacobian from local coordinates to rho/z
  Acts::ActsMatrix<2, 2> jac =
      jacXyzToRhoZ *
      rotLocalToGlobal.block<3, 2>(Acts::ePos0, Acts::ePos0);
  // compute rho/z variance
  Acts::ActsVector<2> var = (jac * localCov * jac.transpose()).diagonal();

  /*
   * From Acts v17 to v19 the scattering uncertainty value allowed was changed.
   * This led to a decrease in efficiency. To offset this, we scale the
   * uncertainties by a tuned factor that gives the v17 performance
   * Track reconstruction is an art as much as it is a science...
   */
  SpacePointPtr spPtr(new SpacePoint{key, x, y, z, r, surf->geometryId(), var[0] * m_uncfactor, var[1] * m_uncfactor});

  if (Verbosity() > 2)
  {
    std::cout << "Space point has "
              << x << ", " << y << ", " << z << " with local coords "
              << localPos.transpose()
              << " with rphi/z variances " << localCov(0, 0)
              << ", " << localCov(1, 1) << " and rotated variances "
              << var[0] << ", " << var[1]
              << " and cluster key "
              << key << " and geo id "
              << surf->geometryId() << std::endl;
  }

  return spPtr;
}

std::vector<const SpacePoint*> PHActsSiliconSeeding::getSiliconSpacePoints(Acts::Extent& rRangeSPExtent,
                                                                           const int strobe)
{
  std::vector<const SpacePoint*> spVec;
  unsigned int numSiliconHits = 0;
  unsigned int totNumSiliconHits = 0;
  std::vector<TrkrDefs::TrkrId> dets = {TrkrDefs::TrkrId::mvtxId};
  if (m_searchInIntt)
  {
    dets.push_back(TrkrDefs::TrkrId::inttId);
  }
  for (const auto& det : dets)
  {
    for (const auto& hitsetkey : m_clusterMap->getHitSetKeys(det))
    {
      if (det == TrkrDefs::TrkrId::mvtxId)
      {
        auto strobeId = MvtxDefs::getStrobeId(hitsetkey);
        if (strobeId != strobe)
        {
          continue;
        }
      }
      auto range = m_clusterMap->getClusters(hitsetkey);
      for (auto clusIter = range.first; clusIter != range.second; ++clusIter)
      {
        const auto cluskey = clusIter->first;
        totNumSiliconHits++;
        if (_iteration_map != nullptr && m_nIteration > 0)
        {
          if (_iteration_map->getIteration(cluskey) < m_nIteration)
          {
            continue;  // skip hits used in a previous iteration
          }
        }

        const auto cluster = clusIter->second;
        const auto hitsetkey_A = TrkrDefs::getHitSetKeyFromClusKey(cluskey);
        const auto surface = m_tGeometry->maps().getSiliconSurface(hitsetkey_A);
        if (!surface)
        {
          continue;
        }

        auto sp = makeSpacePoint(surface, cluskey, cluster).release();
        spVec.push_back(sp);
        rRangeSPExtent.extend({sp->x(), sp->y(), sp->z()});
        numSiliconHits++;
      }
    }
  }
  if (m_seedAnalysis)
  {
    h_nInputMvtxMeas->Fill(numSiliconHits);
  }

  if (Verbosity() > 1)
  {
    std::cout << "Total number of silicon hits : " << totNumSiliconHits << std::endl;
    std::cout << "Seed finding with "
              << numSiliconHits << " hits " << std::endl;
  }
  return spVec;
}

void PHActsSiliconSeeding::configureSPGrid()
{
  m_gridCfg.minPt = m_minSeedPt / m_gridFactor;
  m_gridCfg.rMax = m_rMax;
  m_gridCfg.zMax = m_zMax;
  m_gridCfg.zMin = m_zMin;
  m_gridCfg.deltaRMax = m_deltaRMax;
  m_gridCfg.cotThetaMax = m_cotThetaMax;
  m_gridCfg.impactMax = m_impactMax;
  m_gridCfg.phiBinDeflectionCoverage = m_numPhiNeighbors;

  m_gridOptions.bFieldInZ = m_bField;
  m_gridCfg = m_gridCfg.toInternalUnits();
  m_gridOptions = m_gridOptions.toInternalUnits();
}

Acts::SeedFilterConfig PHActsSiliconSeeding::configureSeedFilter()
{
  Acts::SeedFilterConfig config;
  config.maxSeedsPerSpM = m_maxSeedsPerSpM;
  return config;
}

void PHActsSiliconSeeding::configureSeeder()
{
  /// these are default values that used to be set in Acts
  m_seedFinderCfg.deltaRMinTopSP = 5 * Acts::UnitConstants::mm;
  m_seedFinderCfg.deltaRMaxTopSP = 270 * Acts::UnitConstants::mm;
  m_seedFinderCfg.deltaRMinBottomSP = 5 * Acts::UnitConstants::mm;
  m_seedFinderCfg.deltaRMaxBottomSP = 270 * Acts::UnitConstants::mm;

  /// Limiting location of measurements (e.g. detector constraints)
  m_seedFinderCfg.rMax = m_rMax;
  m_seedFinderCfg.rMin = m_rMin;
  m_seedFinderCfg.zMin = m_zMin;
  m_seedFinderCfg.zMax = m_zMax;

  /// Min/max distance between two measurements in one seed
  m_seedFinderCfg.deltaRMin = m_deltaRMin;
  m_seedFinderCfg.deltaRMax = m_deltaRMax;

  /// Limiting collision region in z
  m_seedFinderCfg.collisionRegionMin = -300. * Acts::UnitConstants::mm;
  m_seedFinderCfg.collisionRegionMax = 300. * Acts::UnitConstants::mm;
  m_seedFinderCfg.sigmaScattering = m_sigmaScattering;
  m_seedFinderCfg.maxSeedsPerSpM = m_maxSeedsPerSpM;
  m_seedFinderCfg.cotThetaMax = m_cotThetaMax;
  m_seedFinderCfg.minPt = m_minSeedPt;
  m_seedFinderOptions.bFieldInZ = m_bField;

  /// Average radiation length traversed per seed
  m_seedFinderCfg.radLengthPerSeed = 0.05;

  /// Maximum impact parameter must be smaller than rMin
  m_seedFinderCfg.impactMax = m_impactMax;

  m_seedFinderCfg.rMinMiddle = 2.5 * Acts::UnitConstants::cm;
  m_seedFinderCfg.rMaxMiddle = 4.3 * Acts::UnitConstants::cm;

  /// Configurations for dealing with misalignment
  m_seedFinderCfg.zAlign = m_zalign;
  m_seedFinderCfg.rAlign = m_ralign;
  m_seedFinderCfg.toleranceParam = m_tolerance;
  m_seedFinderCfg.maxPtScattering = m_maxPtScattering;
  m_seedFinderCfg.sigmaError = m_sigmaError;
  m_seedFinderCfg.helixCut = m_helixcut;

  m_seedFinderCfg =
      m_seedFinderCfg.toInternalUnits().calculateDerivedQuantities();
  m_seedFinderOptions = m_seedFinderOptions.toInternalUnits().calculateDerivedQuantities(m_seedFinderCfg);
}

int PHActsSiliconSeeding::getNodes(PHCompositeNode* topNode)
{
  _cluster_crossing_map = findNode::getClass<TrkrClusterCrossingAssoc>(topNode, "TRKR_CLUSTERCROSSINGASSOC");
  if (!_cluster_crossing_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTERCROSSINGASSOC " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_geomContainerIntt = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  if (!m_geomContainerIntt)
  {
    std::cout << PHWHERE << "CYLINDERGEOM_INTT node not found on node tree"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE << "No ActsGeometry on node tree. Bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (m_useTruthClusters)
  {
    m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode,
                                                            "TRKR_CLUSTER_TRUTH");
  }
  else
  {
    m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode,
                                                            "TRKR_CLUSTER");
  }

  if (!m_clusterMap)
  {
    std::cout << PHWHERE << "No cluster container on the node tree. Bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
int PHActsSiliconSeeding::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHActsSiliconSeeding::createNodes");
  }

  PHNodeIterator dstIter(dstNode);
  PHCompositeNode* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_seedContainer = findNode::getClass<TrackSeedContainer>(topNode, _track_map_name);
  if (!m_seedContainer)
  {
    m_seedContainer = new TrackSeedContainer_v1;
    PHIODataNode<PHObject>* trackNode =
        new PHIODataNode<PHObject>(m_seedContainer, _track_map_name, "PHObject");
    svtxNode->addNode(trackNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsSiliconSeeding::writeHistograms()
{
  m_file->cd();
  h_nInttProj->Write();
  h_nMatchedClusters->Write();
  h_nMvtxHits->Write();
  h_nSeeds->Write();
  h_nActsSeeds->Write();
  h_nTotSeeds->Write();
  h_nInttHits->Write();
  h_nInputMeas->Write();
  h_nHits->Write();
  h_nInputMvtxMeas->Write();
  h_nInputInttMeas->Write();
  h_hits->Write();
  h_zhits->Write();
  h_projHits->Write();
  h_zprojHits->Write();
  h_resids->Write();
  m_file->Write();
  m_file->Close();
}
void PHActsSiliconSeeding::clearTreeVariables()
{
  m_mvtxgx.clear();
  m_mvtxgy.clear();
  m_mvtxgr.clear();
  m_mvtxgz.clear();
  m_projgx = NAN;
  m_projgy = NAN;
  m_projgr = NAN;
  m_projgz = NAN;
  m_projlx = NAN;
  m_projlz = NAN;
  m_clusgx = NAN;
  m_clusgy = NAN;
  m_clusgr = NAN;
  m_clusgz = NAN;
  m_cluslx = NAN;
  m_cluslz = NAN;
}
void PHActsSiliconSeeding::createHistograms()
{
  m_tree = new TTree("seeds", "seed tree");
  m_tree->Branch("seedid", &m_seedid, "m_seedid/I");
  m_tree->Branch("event", &m_event, "m_event/I");
  m_tree->Branch("projgx", &m_projgx, "m_projgx/F");
  m_tree->Branch("projgy", &m_projgy, "m_projgy/F");
  m_tree->Branch("projgr", &m_projgr, "m_projgr/F");
  m_tree->Branch("projgz", &m_projgz, "m_projgz/F");
  m_tree->Branch("projlx", &m_projlx, "m_projlx/F");
  m_tree->Branch("projlz", &m_projlz, "m_projlz/F");
  m_tree->Branch("clusgx", &m_clusgx, "m_clusgx/F");
  m_tree->Branch("clusgy", &m_clusgy, "m_clusgy/F");
  m_tree->Branch("clusgz", &m_clusgz, "m_clusgz/F");
  m_tree->Branch("clusgr", &m_clusgr, "m_clusgr/F");
  m_tree->Branch("cluslx", &m_cluslx, "m_cluslx/F");
  m_tree->Branch("cluslz", &m_cluslz, "m_cluslz/F");
  m_tree->Branch("mvtxgx", &m_mvtxgx);
  m_tree->Branch("mvtxgy", &m_mvtxgy);
  m_tree->Branch("mvtxgz", &m_mvtxgz);
  m_tree->Branch("mvtxgr", &m_mvtxgr);

  h_nMatchedClusters = new TH1F("nMatchedClusters", ";N_{matches}", 50, 0, 50);
  h_nInttProj = new TH2F("nInttProj", ";l_{0}^{proj}-l_{0}^{clus} [cm]; l_{1}^{proj}-l_{1}^{clus} [cm]",
                         10000, -10, 10, 10000, -50, 50);
  h_nMvtxHits = new TH1I("nMvtxHits", ";N_{MVTX}", 6, 0, 6);
  h_nInttHits = new TH1I("nInttHits", ";N_{INTT}", 80, 0, 80);
  h_nHits = new TH2I("nHits", ";N_{MVTX};N_{INTT}", 10, 0, 10,
                     80, 0, 80);
  h_nActsSeeds = new TH1I("nActsSeeds", ";N_{Seeds}", 400, 0, 400);
  h_nSeeds = new TH1I("nActsGoodSeeds", ";N_{Seeds}", 400, 0, 400);
  h_nTotSeeds = new TH1I("nTotSeeds", ";N_{Seeds}", 500, 0, 500);
  h_nInputMeas = new TH1I("nInputMeas", ";N_{Meas}", 2000, 0, 2000);
  h_nInputMvtxMeas = new TH1I("nInputMvtxMeas", ";N_{meas}^{mvtx}",
                              150, 0, 150);
  h_nInputInttMeas = new TH1I("nInputInttMeas", ";N_{meas}^{intt}",
                              150, 0, 150);
  h_hits = new TH2F("hits", ";x [cm]; y [cm]", 1000, -20, 20,
                    1000, -20, 20);
  h_zhits = new TH2F("zhits", ";z [cm]; r [cm]", 1000, -30, 30,
                     1000, -30, 30);
  h_projHits = new TH2F("projhits", ";x [cm]; y [cm]", 1000, -20, 20,
                        1000, -20, 20);
  h_zprojHits = new TH2F("zprojhits", ";z [cm]; r [cm]", 1000, -30, 30,
                         1000, -30, 30);
  h_resids = new TH2F("resids", ";z_{resid} [cm]; rphi_{resid} [cm]",
                      100, -1, 1, 100, -1, 1);
}

double PHActsSiliconSeeding::normPhi2Pi(const double phi)
{
  double returnPhi = phi;
  if (returnPhi < -M_PI)
  {
    returnPhi += 2 * M_PI;
  }
  if (returnPhi > M_PI)
  {
    returnPhi -= 2 * M_PI;
  }
  return returnPhi;
}

void PHActsSiliconSeeding::largeGridSpacing(const bool spacing)
{
  if (!spacing)
  {
    m_gridFactor = 1.;
    m_rMax = 50.;
    m_cotThetaMax = 1.335647;
    m_maxSeedPCA = 0.1;
  }
}

void PHActsSiliconSeeding::printSeedConfigs(Acts::SeedFilterConfig& sfconfig)
{
  // helper function for when updating acts to ensure all seeding parameters
  // are actually the same as before

  std::cout << " --------------- SeedFilterConfig ------------------ " << std::endl;
  std::cout << "deltaInvHelixDiameter = "
            << sfconfig.deltaInvHelixDiameter << std::endl
            << "impactWeightFactor = " << sfconfig.impactWeightFactor
            << std::endl
            << "zOriginWeightFactor = "
            << sfconfig.zOriginWeightFactor << std::endl
            << "compatSeedWeight = " << sfconfig.compatSeedWeight
            << std::endl
            << "deltaRMin = " << sfconfig.deltaRMin
            << std::endl
            << "maxSeedsPerSpM = "
            << sfconfig.maxSeedsPerSpM << std::endl
            << "compatSeedLimit = " << sfconfig.compatSeedLimit
            << std::endl
            << "seedWeightIncrement = "
            << sfconfig.seedWeightIncrement << std::endl
            << "numSeedIncrement = " << sfconfig.numSeedIncrement
            << std::endl;

  std::cout << " --------------- SeedFinderConfig ------------------ " << std::endl;

  std::cout << "deltaRMinTopSp = " << m_seedFinderCfg.deltaRMinTopSP
            << std::endl
            << "deltaRMaxTopSP = " << m_seedFinderCfg.deltaRMaxTopSP
            << std::endl
            << "deltaRMinBottomSP = " << m_seedFinderCfg.deltaRMinBottomSP
            << std::endl
            << "deltaRMaxBottomSP = " << m_seedFinderCfg.deltaRMaxBottomSP
            << std::endl
            << "minpt = " << m_seedFinderCfg.minPt
            << std::endl
            << "cotThetaMax = " << m_seedFinderCfg.cotThetaMax
            << std::endl
            << "deltaRMin = " << m_seedFinderCfg.deltaRMin
            << std::endl
            << "deltaRMax = " << m_seedFinderCfg.deltaRMax
            << std::endl
            << "binSizeR = " << m_seedFinderCfg.binSizeR
            << std::endl
            << "deltaRMiddleMinSPRange = " << m_seedFinderCfg.deltaRMiddleMinSPRange
            << std::endl
            << "deltaRMiddleMaxSPRange = " << m_seedFinderCfg.deltaRMiddleMaxSPRange
            << std::endl
            << "rMinMiddle = " << m_seedFinderCfg.rMinMiddle
            << std::endl
            << "rMaxMiddle = " << m_seedFinderCfg.rMaxMiddle
            << std::endl
            << "deltaZMax = " << m_seedFinderCfg.deltaZMax
            << std::endl
            << "rMax = " << m_seedFinderCfg.rMax
            << std::endl
            << "rMin = " << m_seedFinderCfg.rMin
            << std::endl
            << "zMin = " << m_seedFinderCfg.zMin
            << std::endl
            << "zMax = " << m_seedFinderCfg.zMax
            << std::endl
            << "collisionRegionMin = " << m_seedFinderCfg.collisionRegionMin
            << std::endl
            << "collisionRegionMax = " << m_seedFinderCfg.collisionRegionMax
            << std::endl
            << "sigmaScattering = " << m_seedFinderCfg.sigmaScattering
            << std::endl
            << "maxSeedsPerSpM = " << m_seedFinderCfg.maxSeedsPerSpM
            << std::endl
            << "bFieldInZ = " << m_seedFinderOptions.bFieldInZ
            << std::endl
            << "radLengthPerSeed = " << m_seedFinderCfg.radLengthPerSeed
            << std::endl
            << "impactMax = " << m_seedFinderCfg.impactMax
            << std::endl
            << "zAlign = " << m_seedFinderCfg.zAlign
            << std::endl
            << "rAlign = " << m_seedFinderCfg.rAlign
            << std::endl
            << "toleranceParam = " << m_seedFinderCfg.toleranceParam
            << std::endl
            << "maxPtScattering = " << m_seedFinderCfg.maxPtScattering
            << std::endl
            << "sigmaError = " << m_seedFinderCfg.sigmaError
            << std::endl
            << "helixcut = " << m_seedFinderCfg.helixcut
            << std::endl;

  std::cout << " --------------- SeedFinderOptions ------------------ " << std::endl;
  std::cout << "beamPos = " << m_seedFinderOptions.beamPos.transpose()
            << std::endl
            << "bFieldInZ = " << m_seedFinderOptions.bFieldInZ
            << std::endl
            << "pTPerHelixRadius = " << m_seedFinderOptions.pTPerHelixRadius
            << std::endl
            << "minHelixDiameter2 = " << m_seedFinderOptions.minHelixDiameter2
            << std::endl
            << "pT2perRadius = " << m_seedFinderOptions.pT2perRadius
            << std::endl
            << "sigmapT2perRadius = " << m_seedFinderOptions.sigmapT2perRadius
            << std::endl;

  std::cout << " --------------- SpacePointGridConfig ------------------ " << std::endl;

  std::cout << "minpt = " << m_gridCfg.minPt << std::endl
            << "rMax = " << m_gridCfg.rMax << std::endl
            << "zMax = " << m_gridCfg.zMax << std::endl
            << "zMin = " << m_gridCfg.zMin << std::endl
            << "deltaRMax = " << m_gridCfg.deltaRMax << std::endl
            << "cotThetaMax = " << m_gridCfg.cotThetaMax << std::endl
            << "impactMax = " << m_gridCfg.impactMax << std::endl
            << "phiBinDeflectionCoverage = " << m_gridCfg.phiBinDeflectionCoverage << std::endl
            << "bFieldInZ = " << m_gridOptions.bFieldInZ << std::endl;
}
