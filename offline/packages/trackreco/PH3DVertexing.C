#include "AssocInfoContainer.h"
#include "PH3DVertexing.h"

#include <trackbase_historic/SvtxClusterMap.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertexMap_v1.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

using namespace std;

PH3DVertexing::PH3DVertexing(const std::string& name)
  : SubsysReco(name)
  , _cluster_map(nullptr)
  , _vertex_map(nullptr)
  , _track_map(nullptr)
  , _assoc_container(nullptr)
{
}

int PH3DVertexing::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PH3DVertexing::InitRun(PHCompositeNode* topNode)
{
  return Setup(topNode);
}

int PH3DVertexing::process_event(PHCompositeNode* topNode)
{
  return Process();
}

int PH3DVertexing::End(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PH3DVertexing::Setup(PHCompositeNode* topNode)
{
  int ret = Fun4AllReturnCodes::ABORTRUN;

  ret = CreateNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PH3DVertexing::CreateNodes(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PH3DVertexing::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  _cluster_map = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
  if (!_cluster_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node SvtxClusterMap" << endl;
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
