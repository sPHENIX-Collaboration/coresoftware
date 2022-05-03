#include "PHTrackFitting.h"

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrDefs.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                // for SubsysReco

#include <phool/getClass.h>
#include <phool/phool.h>                       // for PHWHERE

#include <iostream>                            // for operator<<, basic_ostream

using namespace std;

PHTrackFitting::PHTrackFitting(const std::string& name)
  : SubsysReco(name)
  , _cluster_map(nullptr)
  , _hitsets(nullptr)
  , _vertex_map(nullptr)
  , _track_map(nullptr)
  , _track_map_name("SvtxTrackMap")
{
}

int PHTrackFitting::Init(PHCompositeNode* /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackFitting::InitRun(PHCompositeNode* topNode)
{
  return Setup(topNode);
}

int PHTrackFitting::process_event(PHCompositeNode* /*topNode*/)
{
  return Process();
}

int PHTrackFitting::Setup(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackFitting::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  //_cluster_map = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _hitsets = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if(!_hitsets)
    {
      std::cout << PHWHERE << "No hitset container on node tree. Bailing."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  _vertex_map = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!_vertex_map)
  {
    cout << PHWHERE << " ERROR: Can't find SvtxVertexMap." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name);
  if (!_track_map)
  {
    cout << PHWHERE << " ERROR: Can't find SvtxTrackMap." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
