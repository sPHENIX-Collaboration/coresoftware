#include "PH3DVertexing.h"

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>

//#include <trackbase/TrkrClusterContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                // for SubsysReco

#include <phool/getClass.h>
#include <phool/phool.h>                       // for PHWHERE

#include <iostream>                            // for operator<<, basic_ostream

using namespace std;

PH3DVertexing::PH3DVertexing(const std::string& name)
  : SubsysReco(name)
    //  , _cluster_map(nullptr)
  , _vertex_map(nullptr)
  , _track_map(nullptr)
{
}

int PH3DVertexing::InitRun(PHCompositeNode* topNode)
{
  return Setup(topNode);
}

int PH3DVertexing::process_event(PHCompositeNode* /*topNode*/)
{
  return Process();
}

int PH3DVertexing::Setup(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PH3DVertexing::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  //  _cluster_map = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
  /*_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTERS");
  if (!_cluster_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTERS" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  */
  _vertex_map = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!_vertex_map)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_track_map)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
