#include "SortTracksByCrossing.h"

/// Tracking includes

#include <trackbase/TrkrCluster.h>            // for TrkrCluster
#include <trackbase/TrkrDefs.h>               // for cluskey, getLayer, TrkrId
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTrack::C...
#include <trackbase_historic/SvtxTrackMap.h>

#include <globalvertex/SvtxVertex_v2.h>
#include <globalvertex/SvtxVertexMap_v1.h>

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

  
  for(const auto& [trackkey, track] : *_track_map)
    {
      auto crossing = track->get_crossing();
      std::cout << "trackkey " << trackkey << " crossing " << crossing << std::endl;
 
      _track_vertex_crossing_map->addTrackAssoc(crossing, trackkey);
    }

  for(const auto& [vtxkey, vertex] : *_svtx_vertex_map)
    {
       std::cout << "Vertex ID: " << vtxkey << " vertex crossing " << vertex->get_beam_crossing() << " list of tracks: " << std::endl;

       std::set<short int> crossings;
       short int crossing = -1000;
       for (auto trackiter = vertex->begin_tracks(); trackiter != vertex->end_tracks(); ++trackiter)
       {
          SvtxTrack* track = _track_map->get(*trackiter);
          if (!track)
          {
            continue;
          }

          auto siseed = track->get_silicon_seed();
          short int intt_crossing = siseed->get_crossing();

          crossing = track->get_crossing();
          std::cout << " vtxid " << vtxkey  << " crossing " << crossing << " intt_crossing " << intt_crossing 
	  // << " siid " << siid
          << " trackID " << *trackiter
          << " track Z " << track->get_z()
	  << " X " << track->get_x()
	  << " Y " << track->get_y()
	  << " quality " << track->get_quality() 
	  << " pt " << track->get_pt()
          << std::endl;
          siseed->identify();

          crossings.insert(crossing);
       }

      if(crossings.size() > 1)
	{
           std::cout << "Warning: more than one crossing value for vertex ID " << vtxkey << std::endl;
        }

      _track_vertex_crossing_map->addVertexAssoc(crossing, vtxkey);
    }

  // print the results
  _track_vertex_crossing_map->identify();

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

  _svtx_vertex_map = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");
  if (!_svtx_vertex_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SvtxVertexMap: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
 

  return Fun4AllReturnCodes::EVENT_OK;
}
