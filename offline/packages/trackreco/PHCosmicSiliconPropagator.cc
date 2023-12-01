
#include "PHCosmicSiliconPropagator.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrClusterCrossingAssoc.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackSeed_v1.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v1.h>

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
PHCosmicSiliconPropagator::~PHCosmicSiliconPropagator()
{
}

//____________________________________________________________________________..
int PHCosmicSiliconPropagator::Init(PHCompositeNode*)
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
int PHCosmicSiliconPropagator::process_event(PHCompositeNode*)
{
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
        rzpoints.push_back(std::make_pair(r(globPos.x(), globPos.y()), globPos.z()));
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

      for (auto& key : ckeys)
      {
        auto cluster = _cluster_map->findCluster(key);
        auto clusglob = _tgeometry->getGlobalPosition(key, cluster);
        auto pca = TrackFitUtils::get_helix_pca(fitparams, clusglob);
        float dcaz = (pca - clusglob).z();

        if (fabs(dcaz) < _dca_z_cut)
        {
          newClusKeys.push_back(key);
        }
      }
    }

    if (nClusters > 0)
    {
      std::unique_ptr<TrackSeed_v1> si_seed = std::make_unique<TrackSeed_v1>();
      for (auto& key : newClusKeys)
      {
        bool isTpcKey = false;

        if (TrkrDefs::getTrkrId(key) == TrkrDefs::TrkrId::tpcId ||
            TrkrDefs::getTrkrId(key) == TrkrDefs::TrkrId::micromegasId)
        {
          isTpcKey = true;
        }

        if (!isTpcKey)
        {
          si_seed->insert_cluster_key(key);
        }
        else if (isTpcKey)
        {
          tpcseed->insert_cluster_key(key);
        }
      }

      si_seed->circleFitByTaubin(_cluster_map, _tgeometry, 0, 8);
      si_seed->lineFit(_cluster_map, _tgeometry, 0, 8);
      tpcseed->circleFitByTaubin(_cluster_map, _tgeometry, 0, 57);
      tpcseed->lineFit(_cluster_map, _tgeometry, 0, 57);
      TrackSeed* mapped_seed = _si_seeds->insert(si_seed.get());
      std::unique_ptr<SvtxTrackSeed_v1> full_seed = std::make_unique<SvtxTrackSeed_v1>();
      int tpcind = _tpc_seeds->find(tpcseed);
      int siind = _si_seeds->find(mapped_seed);
      full_seed->set_tpc_seed_index(tpcind);
      full_seed->set_silicon_seed_index(siind);
      _svtx_seeds->insert(full_seed.get());
      if (Verbosity() > 3)
      {
        std::cout << "final seeds" << std::endl;
        si_seed->identify();
        tpcseed->identify();
      }
    }
    else
    {
      // no other clusters found, put TPC-only seed in SvtxTrackSeedContainer
      std::unique_ptr<SvtxTrackSeed_v1> partial_seed = std::make_unique<SvtxTrackSeed_v1>();
      int tpc_seed_index = _tpc_seeds->find(tpcseed);
      partial_seed->set_tpc_seed_index(tpc_seed_index);
      _svtx_seeds->insert(partial_seed.get());
    }
  }

  if (Verbosity() > 2)
  {
    std::cout << "svtx seed map size is " << _svtx_seeds->size() << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHCosmicSiliconPropagator::End(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCosmicSiliconPropagator::createSeedContainer(TrackSeedContainer*& container, const std::string &container_name, PHCompositeNode* topNode)
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
