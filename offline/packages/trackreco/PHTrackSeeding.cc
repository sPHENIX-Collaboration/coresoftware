#include "PHTrackSeeding.h"

#include "AssocInfoContainer.h"

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                        // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                      // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>                         // for PHWHERE

#include <iostream>                              // for operator<<, endl

using namespace std;

PHTrackSeeding::PHTrackSeeding(const std::string& name)
  : SubsysReco(name)
{
}

int PHTrackSeeding::InitRun(PHCompositeNode* topNode)
{
  return Setup(topNode);
}

int PHTrackSeeding::process_event(PHCompositeNode* topNode)
{
  return Process(topNode);
}

int PHTrackSeeding::End(PHCompositeNode* topNode)
{
  return End();
}

int PHTrackSeeding::Setup(PHCompositeNode* topNode)
{
  //cout << PHWHERE << "Entering Setup" << endl;
 
  int ret = CreateNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackSeeding::CreateNodes(PHCompositeNode* topNode)
{
  // create nodes...
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = static_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cerr << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  PHNodeIterator iter_dst(dstNode);

  // Create the SVTX node
  PHCompositeNode* tb_node =
      dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode",
                                                        "SVTX"));
  if (!tb_node)
  {
    tb_node = new PHCompositeNode("SVTX");
    dstNode->addNode(tb_node);
    if (Verbosity() > 0)
      cout << PHWHERE << "SVTX node added" << endl;
  }

  
  _track_map = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name);
  if (!_track_map)
    {
      _track_map = new SvtxTrackMap_v1;
      PHIODataNode<PHObject>* tracks_node = 
	new PHIODataNode<PHObject>(_track_map, _track_map_name, "PHObject");
      tb_node->addNode(tracks_node);
      if (Verbosity() > 0){
	cout << PHWHERE << "Svtx/" <<_track_map_name  << " node added" << endl;
      }
    }
  if(Verbosity() > 0)
    _track_map->identify();

  _assoc_container = new AssocInfoContainer;
  PHIODataNode<PHObject>* assoc_node = new PHIODataNode<PHObject>(
      _assoc_container, "AssocInfoContainer", "PHObject");
  tb_node->addNode(assoc_node);
  if (Verbosity() > 0)
    cout << PHWHERE << "Svtx/AssocInfoContainer node added" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackSeeding::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------
  if(_use_truth_clusters)
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_TRUTH");
  else
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  if (!_cluster_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!_cluster_hit_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTERHITASSOC" << endl;
  }

  _vertex_map = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!_vertex_map)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
 
 _track_map = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name);
  if (!_track_map)
  {
    cerr << PHWHERE << " ERROR: Can't find " << _track_map_name << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _assoc_container = findNode::getClass<AssocInfoContainer>(topNode, "AssocInfoContainer");
  if (!_assoc_container)
  {
    cerr << PHWHERE << " ERROR: Can't find AssocInfoContainer." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
