#include "PHTrackSelector.h"

/// Tracking includes

#include <trackbase/TrkrCluster.h>            // for TrkrCluster
#include <trackbase/TrkrDefs.h>               // for cluskey, getLayer, TrkrId
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTrack::C...
#include <trackbase_historic/SvtxTrackMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>                         // for PHIODataNode
#include <phool/PHNode.h>                               // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                             // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <cmath>                              // for sqrt, fabs, atan2, cos
#include <iostream>                           // for operator<<, basic_ostream
#include <map>                                // for map
#include <set>                                // for _Rb_tree_const_iterator
#include <utility>                            // for pair, make_pair

//____________________________________________________________________________..
PHTrackSelector::PHTrackSelector(const std::string &name)
  : SubsysReco(name)
  , _track_map_name("SvtxTrackMap")
  , _n_iter(1)
  , min_tpc_clusters(40)
  , min_mvtx_hits(2)
  , min_intt_hits(1)
  , max_chi2_ndf(30)
{

}

//____________________________________________________________________________..
PHTrackSelector::~PHTrackSelector()
{

}

//____________________________________________________________________________..
int PHTrackSelector::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return ret;
}

//____________________________________________________________________________..
int PHTrackSelector::process_event(PHCompositeNode */*topNode*/)
{

  if(Verbosity() > 0)
    std::cout << PHWHERE << " track map size " << _track_map->size() << std::endl; 

  if(_n_iter <=1) _iteration_map->Reset();

  std::set<unsigned int> track_delete_list;

  unsigned int good_track = 0;  // for diagnostic output only
  unsigned int bad_track = 0;   // tracks to keep

  // make a list of tracks that did not make the keep list
  for(auto track_it = _track_map->begin(); track_it != _track_map->end(); ++track_it)
    {
      auto id = track_it->first;
      _track = track_it->second;
      bool delete_track = false;
      double chi2_ndf = _track->get_chisq() / _track->get_ndf();
      int ntpc = 0;
      int nintt = 0;
      int nmvtx = 0;

      for (SvtxTrack::ConstClusterKeyIter iter = _track->begin_cluster_keys();
	   iter != _track->end_cluster_keys();
	   ++iter){
	TrkrDefs::cluskey cluster_key = *iter;
	//	unsigned int layer = TrkrDefs::getLayer(cluster_key);

	if(TrkrDefs::getTrkrId(cluster_key) == TrkrDefs::mvtxId ) nmvtx++;
	if(TrkrDefs::getTrkrId(cluster_key) == TrkrDefs::inttId ) nintt++;
	if(TrkrDefs::getTrkrId(cluster_key) == TrkrDefs::tpcId )  ntpc++;

      }

      if(chi2_ndf > max_chi2_ndf) delete_track = true;
      if(ntpc < min_tpc_clusters) delete_track = true;
      if(ntpc < min_tpc_clusters) delete_track = true;
      if(nintt < min_intt_hits) delete_track = true;
      if(nmvtx < min_mvtx_hits) delete_track = true;

      if(delete_track){
	track_delete_list.insert(id);
	bad_track++;
      }
      else{
	//add clusters to cluster used list
	for (SvtxTrack::ConstClusterKeyIter iter = _track->begin_cluster_keys();
	     iter != _track->end_cluster_keys();
	     ++iter){
	  TrkrDefs::cluskey cluster_key = *iter;
	  _iteration_map->addIteration(cluster_key,_n_iter);
	}
	good_track++;
      }
      
    }
  if(Verbosity() > 0)
    std::cout << " Number of good tracks is " << good_track << " bad tracks " << bad_track << std::endl; 
  
  // delete failed tracks
  for(auto it = track_delete_list.begin(); it != track_delete_list.end(); ++it){
    if(Verbosity() > 1)
      std::cout << " erasing track ID " << *it << std::endl;
    _track_map->erase(*it);
  }
  
  if(Verbosity() > 0)
    std::cout << " track_delete_list size " << track_delete_list.size() << std::endl;

  // delete failed tracks
  for(auto it = track_delete_list.begin(); it != track_delete_list.end(); ++it)
    {
      if(Verbosity() > 1)
	std::cout << " erasing track ID " << *it << std::endl;
      _track_map->erase(*it);
    }

  if(Verbosity() > 0)
    std::cout << "Track map size after choosing best silicon match: " << _track_map->size() << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackSelector::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int  PHTrackSelector::GetNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  //  TrackClusterIterationMapv1* _iteration_map;
  _track_map = findNode::getClass<SvtxTrackMap>(topNode,  _track_map_name);
  if (!_track_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Create the Cluster Iteration Map node if required
  _iteration_map= findNode::getClass<TrkrClusterIterationMapv1>(topNode, "CLUSTER_ITERATION_MAP");
  if (!_iteration_map)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode =
      dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    _iteration_map = new TrkrClusterIterationMapv1;
    PHIODataNode<PHObject> *TrkrClusterIterationMapv1 =
      new PHIODataNode<PHObject>(_iteration_map, "CLUSTER_ITERATION_MAP", "PHObject");
    DetNode->addNode(TrkrClusterIterationMapv1);
  }



  return Fun4AllReturnCodes::EVENT_OK;
}
