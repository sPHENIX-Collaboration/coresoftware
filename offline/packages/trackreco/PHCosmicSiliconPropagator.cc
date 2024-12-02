
#include "PHCosmicSiliconPropagator.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterCrossingAssoc.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackSeed_v1.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v2.h>
#include <trackbase_historic/TrackSeedHelper.h>

#include <cmath>

namespace
{
  template <class T>
  inline constexpr T square(T& x)
  {
    return x * x;
  }
  template <class T>
  inline constexpr T r(T& x, T& y)
  {
    return std::sqrt(square(x) + square(y));
  }

}  // namespace
//____________________________________________________________________________..
PHCosmicSiliconPropagator::PHCosmicSiliconPropagator(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
PHCosmicSiliconPropagator::~PHCosmicSiliconPropagator() = default;

//____________________________________________________________________________..
int PHCosmicSiliconPropagator::Init(PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHCosmicSiliconPropagator::InitRun(PHCompositeNode* topNode)
{
  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!_tgeometry)
  {
    std::cout << "No Acts tracking geometry, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _tpc_seeds = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!_tpc_seeds)
  {
    std::cout << "No TpcTrackSeedContainer, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _si_seeds = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if (!_si_seeds)
  {
    std::cout << "No SiliconTrackSeedContainer, creating..." << std::endl;
    if (createSeedContainer(_si_seeds, "SiliconTrackSeedContainer", topNode) != Fun4AllReturnCodes::EVENT_OK)
    {
      std::cout << "Cannot create, exiting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  _svtx_seeds = findNode::getClass<TrackSeedContainer>(topNode, _track_map_name);
  if (!_svtx_seeds)
  {
    std::cout << "No " << _track_map_name << " found, creating..." << std::endl;
    if (createSeedContainer(_svtx_seeds, _track_map_name, topNode) != Fun4AllReturnCodes::EVENT_OK)
    {
      std::cout << "Cannot create, exiting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHCosmicSiliconPropagator::process_event(PHCompositeNode* /*unused*/)
{
  if (m_resetContainer)
  {
    _svtx_seeds->Reset();
  }

  for (auto& tpcseed : *_tpc_seeds)
  {
    if (!tpcseed)
    {
      continue;
    }
    if (Verbosity() > 3)
    {
      std::cout << "Tpc seed to start out is " << std::endl;
      tpcseed->identify();
    }
    std::vector<Acts::Vector3> tpcClusPos;
    std::vector<TrkrDefs::cluskey> tpcClusKeys;

    std::copy(tpcseed->begin_cluster_keys(), tpcseed->end_cluster_keys(),
              std::back_inserter(tpcClusKeys));

    TrackFitUtils::getTrackletClusters(_tgeometry, _cluster_map,
                                       tpcClusPos, tpcClusKeys);
    std::vector<TrkrDefs::cluskey> newClusKeys;
    std::vector<Acts::Vector3> newClusPos;
    unsigned int nClusters = 0;

    if (m_zeroField)
    {
      //! Line fit both in x-y and r-z
      TrackFitUtils::position_vector_t rzpoints, xypoints;
      for (auto& globPos : tpcClusPos)
      {
        xypoints.push_back(std::make_pair(globPos.x(), globPos.y()));
        float clusr = r(globPos.x(), globPos.y());
        if (globPos.y() < 0)
        {
          clusr *= -1;
        }
        rzpoints.push_back(std::make_pair(globPos.z(), clusr));
      }

      auto xyLineParams = TrackFitUtils::line_fit(xypoints);
      auto rzLineParams = TrackFitUtils::line_fit(rzpoints);

      std::vector<TrkrDefs::cluskey> newClusKeysrz;
      std::vector<Acts::Vector3> newClusPosrz;
      std::vector<TrkrDefs::cluskey> newClusKeysxy;
      std::vector<Acts::Vector3> newClusPosxy;
      // now add clusters along lines
      nClusters = TrackFitUtils::addClustersOnLine(xyLineParams,
                                                   true,
                                                   _dca_xy_cut,
                                                   _tgeometry,
                                                   _cluster_map,
                                                   newClusPosxy,
                                                   newClusKeysxy,
                                                   0, 56);
      int nrzClusters = TrackFitUtils::addClustersOnLine(rzLineParams,
                                                         false,
                                                         _dca_z_cut,
                                                         _tgeometry,
                                                         _cluster_map,
                                                         newClusPosrz,
                                                         newClusKeysrz,
                                                         0, 56);

      if (Verbosity() > 3)
      {
        std::cout << "nrz clus " << nrzClusters << " and nxy clusters " << nClusters
                  << std::endl;
      }
      std::set_intersection(newClusKeysxy.begin(), newClusKeysxy.end(),
                            newClusKeysrz.begin(), newClusKeysrz.end(), std::back_inserter(newClusKeys));
      if (m_resetContainer)
      {
        for (auto& keys : {newClusKeysxy, newClusKeysrz})
        {
          for (auto& key : keys)
          {
            if (TrkrDefs::getTrkrId(key) == TrkrDefs::TrkrId::micromegasId)
            {
              newClusKeys.push_back(key);
            }
          }
        }
      }

      if (Verbosity() > 3)
      {
        for (auto key : newClusKeysxy)
        {
          auto cluster = _cluster_map->findCluster(key);
          auto clusglob = _tgeometry->getGlobalPosition(key, cluster);
          std::cout << "Found key " << key << " for xy cosmic in layer " << (unsigned int) TrkrDefs::getLayer(key)
                    << " with pos " << clusglob.transpose() << std::endl;
        }
        for (auto key : newClusKeysrz)
        {
          auto cluster = _cluster_map->findCluster(key);
          auto clusglob = _tgeometry->getGlobalPosition(key, cluster);
          std::cout << "Found key " << key << " for rz cosmic in layer " << (unsigned int) TrkrDefs::getLayer(key)
                    << " with pos " << clusglob.transpose() << std::endl;
        }
      }
    }
    else
    {
      std::vector<float> fitparams = TrackFitUtils::fitClusters(tpcClusPos, tpcClusKeys);
      //! There weren't enough clusters to fit
      if (fitparams.size() == 0)
      {
        continue;
      }

      std::vector<TrkrDefs::cluskey> ckeys;
      nClusters = TrackFitUtils::addClusters(fitparams, _dca_xy_cut, _tgeometry, _cluster_map,
                                             newClusPos, ckeys, 0, 56);
      TrackFitUtils::position_vector_t yzpoints;
      for (auto& globPos : tpcClusPos)
      {
        yzpoints.push_back(std::make_pair(globPos.y(), globPos.z()));
      }

      auto yzLineParams = TrackFitUtils::line_fit(yzpoints);
      float yzslope = std::get<0>(yzLineParams);
      float yzint = std::get<1>(yzLineParams);
      for (auto& key : ckeys)
      {
        auto cluster = _cluster_map->findCluster(key);
        auto clusglob = _tgeometry->getGlobalPosition(key, cluster);

        float projz = clusglob.y() * yzslope + yzint;

        if (std::fabs(projz - clusglob.z()) < _dca_z_cut)
        {
          newClusKeys.push_back(key);
        }
      }
    }
    //! only keep long seeds
    if ((tpcClusKeys.size() + newClusKeys.size() > 25))
    {
      // TODO: should include distortion corrections
      std::unique_ptr<TrackSeed_v2> si_seed = std::make_unique<TrackSeed_v2>();
      std::map<TrkrDefs::cluskey, Acts::Vector3> silposmap, tpcposmap;
      for (auto& key : tpcClusKeys)
      {
        auto cluster = _cluster_map->findCluster(key);
        auto clusglob = _tgeometry->getGlobalPosition(key, cluster);
        tpcposmap.emplace(key, clusglob);
      }
      for (auto& key : newClusKeys)
      {
        bool isTpcKey = false;
        auto cluster = _cluster_map->findCluster(key);
        auto clusglob = _tgeometry->getGlobalPosition(key, cluster);
        if (TrkrDefs::getTrkrId(key) == TrkrDefs::TrkrId::tpcId ||
            TrkrDefs::getTrkrId(key) == TrkrDefs::TrkrId::micromegasId)
        {
          isTpcKey = true;
        }
        if (!isTpcKey)
        {
          si_seed->insert_cluster_key(key);
          silposmap.emplace(key, clusglob);
        }
        else
        {
          tpcseed->insert_cluster_key(key);
          tpcposmap.emplace(key, clusglob);
        }
      }

      TrackSeedHelper::circleFitByTaubin(si_seed.get(), silposmap, 0, 8);
      TrackSeedHelper::lineFit(si_seed.get(), silposmap);

      TrackSeedHelper::circleFitByTaubin(tpcseed, tpcposmap, 7, 57);
      TrackSeedHelper::lineFit(tpcseed, tpcposmap, 7, 57);

      TrackSeed* mapped_seed = _si_seeds->insert(si_seed.get());
      std::unique_ptr<SvtxTrackSeed_v1> full_seed = std::make_unique<SvtxTrackSeed_v1>();
      int tpcind = _tpc_seeds->find(tpcseed);
      int siind = _si_seeds->find(mapped_seed);
      full_seed->set_tpc_seed_index(tpcind);
      if (si_seed->size_cluster_keys() > 0)
      {
        full_seed->set_silicon_seed_index(siind);
      }
      _svtx_seeds->insert(full_seed.get());
      if (Verbosity() > 3)
      {
        std::cout << "final seeds" << std::endl;
        si_seed->identify();
        tpcseed->identify();
      }
    }
  }

  if (Verbosity() > 2)
  {
    std::cout << "svtx seed map size is " << _svtx_seeds->size() << std::endl;
    int i = 0;
    for (auto& seed : *_svtx_seeds)
    {
      std::cout << "seed " << i << " is composed of " << std::endl;
      _tpc_seeds->get(seed->get_tpc_seed_index())->identify();
      if (_si_seeds->get(seed->get_silicon_seed_index()))
      {
        _si_seeds->get(seed->get_silicon_seed_index())->identify();
      }
      ++i;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHCosmicSiliconPropagator::End(PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCosmicSiliconPropagator::createSeedContainer(TrackSeedContainer*& container, const std::string& container_name, PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHSiliconHelicalPropagator::createNodes");
  }

  PHNodeIterator dstIter(dstNode);
  PHCompositeNode* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  container = findNode::getClass<TrackSeedContainer>(topNode, container_name);
  if (!container)
  {
    container = new TrackSeedContainer_v1;
    PHIODataNode<PHObject>* trackNode = new PHIODataNode<PHObject>(container, container_name, "PHObject");
    svtxNode->addNode(trackNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
