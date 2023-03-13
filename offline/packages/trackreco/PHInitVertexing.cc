#include "PHInitVertexing.h"

#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertexMap_v1.h>
#include <trackbase/TrkrClusterContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                   // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                         // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                       // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>                          // for PHWHERE

#include <iostream>                               // for operator<<, endl

using namespace std;

PHInitVertexing::PHInitVertexing(const std::string& name)
  : SubsysReco(name)
{}

int PHInitVertexing::InitRun(PHCompositeNode* topNode)
{
  return Setup(topNode);
}

int PHInitVertexing::process_event(PHCompositeNode* topNode)
{
  return Process(topNode);
}

int PHInitVertexing::Setup(PHCompositeNode* topNode)
{
  int ret = CreateNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHInitVertexing::CreateNodes(PHCompositeNode* topNode)
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

  _vertex_map = new SvtxVertexMap_v1;
  PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(
      _vertex_map, "SvtxVertexMap", "PHObject");
  tb_node->addNode(vertexes_node);
  if (Verbosity() > 0)
    cout << PHWHERE << "Svtx/SvtxVertexMap node added" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHInitVertexing::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node TrkrClusterContainer" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Pull the reconstructed track information off the node tree...
  _vertex_map = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!_vertex_map)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
