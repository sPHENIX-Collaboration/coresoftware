#include "PHSiliconHelicalPropagator.h"

#include <trackbase_historic/TrackSeed_v1.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

PHSiliconHelicalPropagator::PHSiliconHelicalPropagator(std::string name) : SubsysReco(name)
{
  _fitter = new HelicalFitter();
}

PHSiliconHelicalPropagator::~PHSiliconHelicalPropagator()
{
  delete _fitter;
}

int PHSiliconHelicalPropagator::InitRun(PHCompositeNode* topNode)
{
  _fitter->InitRun(topNode);

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if(!_tgeometry)
  {
    std::cout << "No Acts tracking geometry, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _si_seeds = findNode::getClass<TrackSeedContainer>(topNode,"SiliconTrackSeedContainer");
  if(!_si_seeds)
  {
    std::cout << "No SiliconTrackSeedContainer, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _track_map = findNode::getClass<SvtxTrackMap>(topNode,_track_map_name);
  if(!_track_map)
  {
    std::cout << "No track map found, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSiliconHelicalPropagator::process_event(PHCompositeNode* /*topNode*/)
{
  for(SvtxTrackMap::Iter track_iter = _track_map->begin(); track_iter != _track_map->end(); ++track_iter)
  {
    SvtxTrack* track = track_iter->second;
    if(!track) continue;

    std::vector<Acts::Vector3> clusterPositions;
    std::vector<TrkrDefs::cluskey> clusterKeys;
    TrackSeed* tpcseed = track->get_tpc_seed();
    _fitter->getTrackletClusters(tpcseed,clusterPositions,clusterKeys);
    std::vector<float> fitparams = _fitter->fitClusters(clusterPositions,clusterKeys);
    
    std::vector<TrkrDefs::cluskey> si_clusterKeys;
    std::vector<Acts::Vector3> si_clusterPositions;
    unsigned int nSiClusters = _fitter->addSiliconClusters(fitparams,si_clusterPositions,si_clusterKeys);
    if(Verbosity()>0) std::cout << "adding " << nSiClusters << " cluster keys" << std::endl;
    
    std::unique_ptr<TrackSeed_v1> si_seed = std::make_unique<TrackSeed_v1>();
    for(auto clusterkey : si_clusterKeys)
    {
      si_seed->insert_cluster_key(clusterkey);
    }
    si_seed->circleFitByTaubin(_cluster_map,_tgeometry,0,8);
    si_seed->lineFit(_cluster_map,_tgeometry,0,8);
    
    TrackSeed* mapped_seed = _si_seeds->insert(si_seed.get());
    
    track->set_silicon_seed(mapped_seed);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSiliconHelicalPropagator::End(PHCompositeNode* /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
