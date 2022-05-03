#include "PHTrackSetMerging.h"

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <trackbase/TrkrClusterContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                        // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                      // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>                         // for PHWHERE

#include <iostream>                              // for operator<<, basic_os...

using namespace std;

PHTrackSetMerging::PHTrackSetMerging(const std::string& name)
  : SubsysReco(name)
{
}

int PHTrackSetMerging::Init(PHCompositeNode* /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackSetMerging::InitRun(PHCompositeNode* topNode)
{
  return Setup(topNode);
}

int PHTrackSetMerging::process_event(PHCompositeNode* topNode)
{
  return Process(topNode);
}

int PHTrackSetMerging::End(PHCompositeNode* /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;//End();
}

int PHTrackSetMerging::Setup(PHCompositeNode* topNode)
{
  int ret = CreateNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackSetMerging::CreateNodes(PHCompositeNode* topNode)
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
      cout << "SVTX node added" << endl;
  }

  _track_map_out = new SvtxTrackMap_v1;
  PHIODataNode<PHObject>* tracks_node = new PHIODataNode<PHObject>(
      _track_map_out, _track_map_name_out, "PHObject");

  //  PHIODataNode<PHObject>* tracks_node = new PHIODataNode<PHObject>(
  //   _track_map_out, "SvtxTrackMapOut", "PHObject");

  tb_node->addNode(tracks_node);
  if (Verbosity() > 0){
    cout << "Svtx/SvtxTrackMapOut node added" << endl;
    // cout << "Svtx/" << _track_map_name << " node added" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackSetMerging::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _vertex_map = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!_vertex_map)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  
  _track_map_in1 = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name_in1);
  if (!_track_map_in1)
  {
    cerr << PHWHERE << " ERROR: Can't find " << _track_map_name_in1 << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map_in2 = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name_in2);
  if (!_track_map_in2)
  {
    cerr << PHWHERE << " ERROR: Can't find " << _track_map_name_in2 << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map_out = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name_out);
  if (!_track_map_out)
  {
    cerr << PHWHERE << " ERROR: Can't find " << _track_map_name_out << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
