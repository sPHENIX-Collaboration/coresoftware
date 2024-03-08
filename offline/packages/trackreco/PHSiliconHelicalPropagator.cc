#include "PHSiliconHelicalPropagator.h"

#include <trackbase_historic/SvtxTrackSeed_v1.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v1.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

namespace 
{
  template <typename T> int sgn(const T& x)
  { 
    if(x > 0) return 1;
    else return -1;
  }
}

PHSiliconHelicalPropagator::PHSiliconHelicalPropagator(const std::string &name)
  : SubsysReco(name)
{
}

PHSiliconHelicalPropagator::~PHSiliconHelicalPropagator()
{
}

int PHSiliconHelicalPropagator::InitRun(PHCompositeNode* topNode)
{
  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _cluster_crossing_map = findNode::getClass<TrkrClusterCrossingAssoc>(topNode, "TRKR_CLUSTERCROSSINGASSOC");
  if (!_cluster_crossing_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find TRKR_CLUSTERCROSSINGASSOC " << std::endl;
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

int PHSiliconHelicalPropagator::createSeedContainer(TrackSeedContainer*& container, const std::string &container_name, PHCompositeNode* topNode)
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

int PHSiliconHelicalPropagator::process_event(PHCompositeNode* /*topNode*/)
{
  for (unsigned int seedID = 0; seedID < _tpc_seeds->size(); ++seedID)
  {
    TrackSeed* tpcseed = _tpc_seeds->get(seedID);
    if (!tpcseed) continue;

    std::vector<Acts::Vector3> clusterPositions;
    std::vector<TrkrDefs::cluskey> clusterKeys;
    for (auto iter = tpcseed->begin_cluster_keys();
         iter != tpcseed->end_cluster_keys(); ++iter)
    {
      clusterKeys.push_back(*iter);
    }
    TrackFitUtils::getTrackletClusters(_tgeometry, _cluster_map, clusterPositions, clusterKeys);
    std::vector<float> fitparams = TrackFitUtils::fitClusters(clusterPositions, clusterKeys);
    //! There weren't enough clusters to fit
    if (fitparams.size() == 0)
    {
      continue;
    }
    std::vector<TrkrDefs::cluskey> si_clusterKeys;
    std::vector<Acts::Vector3> si_clusterPositions;

    unsigned int nSiClusters = TrackFitUtils::addClusters(fitparams, 1000., _tgeometry, _cluster_map, si_clusterPositions, si_clusterKeys,0,6);

    if (nSiClusters > 0)
    {
      
      std::unique_ptr<TrackSeed_v1> si_seed = std::make_unique<TrackSeed_v1>();
      std::map<short, int> crossing_frequency;

      Acts::Vector3 layer0global;

      for (auto clusterkey : si_clusterKeys)
      {
	if(TrkrDefs::getLayer(clusterkey) == 0)
	  {
	    auto cluster = _cluster_map->findCluster(clusterkey);
	    layer0global = _tgeometry->getGlobalPosition(clusterkey, cluster);
	  }
	if(TrkrDefs::getTrkrId(clusterkey) == TrkrDefs::mvtxId)
	  {
	    si_seed->insert_cluster_key(clusterkey);
	  }
        else if (TrkrDefs::getTrkrId(clusterkey) == TrkrDefs::inttId)
        {
          auto hit_crossings = _cluster_crossing_map->getCrossings(clusterkey);
          for (auto iter = hit_crossings.first; iter != hit_crossings.second; ++iter)
          {
            short crossing = iter->second;
            if (crossing_frequency.count(crossing) == 0)
              crossing_frequency.insert({crossing, 1});
            else
              crossing_frequency[crossing]++;
          }

	  //! Check that the INTT clusters are in the same quadrant
	  auto cluster = _cluster_map->findCluster(clusterkey);
	  auto global = _tgeometry->getGlobalPosition(clusterkey, cluster);
	  if(sgn(global.x()) == sgn(layer0global.x()) && sgn(global.y()) == sgn(layer0global.y()))
	    {
	      si_seed->insert_cluster_key(clusterkey);
	    }
	}

      }

      if (crossing_frequency.size() > 0)
      {
        short most_common_crossing = (std::max_element(crossing_frequency.begin(), crossing_frequency.end(),
                                                       [](auto entry1, auto entry2)
                                                       { return entry1.second > entry2.second; }))
                                         ->first;
        si_seed->set_crossing(most_common_crossing);
      }
      si_seed->circleFitByTaubin(_cluster_map, _tgeometry, 0, 8);
      si_seed->lineFit(_cluster_map, _tgeometry, 0, 8);

      TrackSeed* mapped_seed = _si_seeds->insert(si_seed.get());

      std::unique_ptr<SvtxTrackSeed_v1> full_seed = std::make_unique<SvtxTrackSeed_v1>();
      int tpc_seed_index = _tpc_seeds->find(tpcseed);
      int si_seed_index = _si_seeds->find(mapped_seed);
      if (Verbosity() > 0)
      {
        std::cout << "inserted " << nSiClusters << " silicon clusters for tpc seed " << tpc_seed_index << std::endl;
        std::cout << "new silicon seed index: " << si_seed_index << std::endl;
      }
      full_seed->set_tpc_seed_index(tpc_seed_index);
      full_seed->set_silicon_seed_index(si_seed_index);
      _svtx_seeds->insert(full_seed.get());
    }
    else
    {
      // no Si clusters found, put TPC-only seed in SvtxTrackSeedContainer
      std::unique_ptr<SvtxTrackSeed_v1> partial_seed = std::make_unique<SvtxTrackSeed_v1>();
      int tpc_seed_index = _tpc_seeds->find(tpcseed);
      partial_seed->set_tpc_seed_index(tpc_seed_index);
      _svtx_seeds->insert(partial_seed.get());
    }
  }
  if (Verbosity() > 2)
    std::cout << "svtx seed map size is " << _svtx_seeds->size() << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSiliconHelicalPropagator::End(PHCompositeNode* /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
