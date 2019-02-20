#include "PHTrackPropagating.h"

#include "AssocInfoContainer.h"

//#include <trackbase_historic/SvtxClusterMap.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertexMap_v1.h>

#include <trackbase/TrkrClusterContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

using namespace std;

PHTrackPropagating::PHTrackPropagating(const std::string& name)
  : SubsysReco(name)
  , _cluster_map(nullptr)
  , _vertex_map(nullptr)
  , _track_map(nullptr)
  , _assoc_container(nullptr)
{
}

int PHTrackPropagating::InitRun(PHCompositeNode* topNode)
{
  return Setup(topNode);
}

int PHTrackPropagating::process_event(PHCompositeNode* topNode)
{
  return Process();
}

int PHTrackPropagating::End(PHCompositeNode* topNode)
{
  End();
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackPropagating::Setup(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackPropagating::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------


  //_cluster_map = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
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

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_track_map)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap." << endl;
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
