#include "SortTracksByCrossing.h"

/// Tracking includes

#include <trackbase/TrkrCluster.h>            // for TrkrCluster
#include <trackbase/TrkrDefs.h>               // for cluskey, getLayer, TrkrId
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTrack::C...
#include <trackbase_historic/SvtxTrackMap.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cmath>                              // for sqrt, fabs, atan2, cos
#include <iostream>                           // for operator<<, basic_ostream
#include <map>                                // for map
#include <set>                                // for _Rb_tree_const_iterator
#include <utility>                            // for pair, make_pair
#include <iomanip>

#include <vector>
#include <cassert>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <functional>

#include <Eigen/Dense>

//____________________________________________________________________________..
SortTracksByCrossing::SortTracksByCrossing(const std::string &name)
  : SubsysReco(name)
{

}

//____________________________________________________________________________..
SortTracksByCrossing::~SortTracksByCrossing()
{

}

//____________________________________________________________________________..
int SortTracksByCrossing::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  ret = CreateNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  return ret;
}

//____________________________________________________________________________..
int SortTracksByCrossing::process_event(PHCompositeNode */*topNode*/)
{

  if(Verbosity() > 0)
    std::cout << PHWHERE << " track map size " << _track_map->size()  << std::endl;

  std::set<short int> crossings;
  
  for(const auto& [trackkey, track] : *_track_map)
    {
      auto crossing = track->get_crossing();
      crossings.insert(crossing);
 
      _track_vertex_crossing_map->addTrackAssoc(crossing, trackkey);

       if(Verbosity() > 0) { std::cout << "trackkey " << trackkey << " crossing " << crossing << std::endl; }
    }

  return Fun4AllReturnCodes::EVENT_OK;
}



int SortTracksByCrossing::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int SortTracksByCrossing::CreateNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  /// Check that it is there
  if (!dstNode)
  {
    std::cerr << "DST Node missing, quitting" << std::endl;
    throw std::runtime_error("failed to find DST node in PHActsInitialVertexFinder::createNodes");
  }

  /// Get the tracking subnode
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

_track_vertex_crossing_map = findNode::getClass<TrackVertexCrossingAssoc>(topNode,"TrackVertexCrossingAssocMap");
  
  if(!_track_vertex_crossing_map)
    {
      _track_vertex_crossing_map = new TrackVertexCrossingAssoc_v1;
      PHIODataNode<PHObject>* trackvertexNode = new PHIODataNode<PHObject>( 
		   _track_vertex_crossing_map, "TrackVertexCrossingAssocMap","PHObject");

      svtxNode->addNode(trackvertexNode);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}
int SortTracksByCrossing::GetNodes(PHCompositeNode* topNode)
{

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_track_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }


  return Fun4AllReturnCodes::EVENT_OK;
}
