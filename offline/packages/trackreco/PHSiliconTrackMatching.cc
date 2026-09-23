#include "PHSiliconTrackMatching.h"

/// Tracking includes
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterCrossingAssoc.h>
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrDefs.h>  // for cluskey, getTrkrId, tpcId

#include <trackbase_historic/SvtxTrackSeed_v2.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v2.h>
#include <trackbase_historic/TrackSeedHelper.h>

#include <globalvertex/SvtxVertex.h>  // for SvtxVertex
#include <globalvertex/SvtxVertexMap.h>

#include <g4main/PHG4Hit.h>       // for PHG4Hit
#include <g4main/PHG4HitDefs.h>   // for keytype
#include <g4main/PHG4Particle.h>  // for PHG4Particle

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/sphenix_constants.h>

#include <TF1.h>
#include <TFile.h>
#include <TNtuple.h>

#include <climits>   // for UINT_MAX
#include <cmath>     // for fabs, sqrt
#include <iostream>  // for operator<<, basic_ostream
#include <memory>
#include <set>      // for _Rb_tree_const_iterator
#include <utility>  // for pair

using namespace std;

//____________________________________________________________________________..
PHSiliconTrackMatching::PHSiliconTrackMatching(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
PHSiliconTrackMatching::~PHSiliconTrackMatching() = default;

//____________________________________________________________________________..
int PHSiliconTrackMatching::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }
  return ret;
}


//____________________________________________________________________________..
int PHSiliconTrackMatching::process_event(PHCompositeNode * /*unused*/)
{
  if(Verbosity() > 2)
  {
    std::cout << " Warning: PHSiliconTpcTrackMatching "
      << ( _zero_field ? "zero field is ON" : " zero field is OFF") << std::endl;
  }
  // _track_map contains the TPC seed track stubs
  // _track_map_silicon contains the silicon seed track stubs
  // _svtx_seed_map contains the combined silicon and tpc track seeds

    // in case these objects are in the input file, we clear the nodes and replace them
  _svtx_seed_map->Reset(); 
  _track_map->Reset();
  
  if (Verbosity() > 0)
  {
    cout << PHWHERE << " TPC track map size " << _track_map->size() << " Silicon track map size " << _track_map_silicon->size() << endl;
  }

  if (_track_map_silicon->empty())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  unsigned int si_id = 0;
  // loop over the silicon seeds and add the crossing to them
  for (unsigned int trackid = 0; trackid != _track_map_silicon->size(); ++trackid)
  {
    _tracklet_si = _track_map_silicon->get(trackid);
    if (!_tracklet_si)
    {
      continue;
    }
    auto crossing = _tracklet_si->get_crossing();
    if (Verbosity() > 8)
    {
      std::cout << " silicon stub: " << trackid << " eta " << _tracklet_si->get_eta()
        << " pt " << _tracklet_si->get_pt() << " si z " << TrackSeedHelper::get_z(_tracklet_si)
        << " crossing " << crossing << std::endl;
    }

    if (Verbosity() > 1)
    {
      cout << " Si track " << trackid << " crossing " << crossing << endl;
    }
    auto dummy = std::make_unique<TrackSeed_v2>();
    dummy->set_qOverR(_tracklet_si->get_qOverR());
    dummy->set_phi(_tracklet_si->get_phi());
    dummy->set_X0(_tracklet_si->get_X0());
    dummy->set_Y0(_tracklet_si->get_Y0());
    dummy->set_Z0(_tracklet_si->get_Z0());
    dummy->set_slope(_tracklet_si->get_slope());
    
    auto svtxseed = std::make_unique<SvtxTrackSeed_v2>();
    svtxseed->set_silicon_seed_index(si_id);
    svtxseed->set_tpc_seed_index(si_id);
    // In pp mode, if a matched track does not have INTT clusters we have to find the crossing geometrically
    // Record the geometrically estimated crossing in the track seeds for later use if needed
    svtxseed->set_crossing_estimate(crossing);
    _track_map->insert(dummy.get());
    _svtx_seed_map->insert(svtxseed.get());

    si_id++;
    if (Verbosity() > 1)
    {
      std::cout << "  combined seed id " << _svtx_seed_map->size() - 1 << " si id " << si_id << " tpc id " << si_id << " crossing estimate " << crossing << std::endl;
    }
  }


  if (Verbosity() > 0)
  {
    std::cout << "final svtx seed map size " << _svtx_seed_map->size() << std::endl;
  }

  if (Verbosity() > 1)
  {
    for (const auto &seed : *_svtx_seed_map)
    {
      seed->identify();
      std::cout << std::endl;
    }

    cout << "PHSiliconTrackMatching::process_event(PHCompositeNode *topNode) Leaving process_event" << endl;
  }
  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}


int PHSiliconTrackMatching::End(PHCompositeNode * /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSiliconTrackMatching::GetNodes(PHCompositeNode *topNode)
{
  //---------------------------------
  // Get additional objects off the Node Tree
  //---------------------------------

  _track_map_silicon = findNode::getClass<TrackSeedContainer>(topNode, _silicon_track_map_name);
  if (!_track_map_silicon)
  {
    cerr << PHWHERE << " ERROR: Can't find SiliconTrackSeedContainer " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map = findNode::getClass<TrackSeedContainer>(topNode, _track_map_name);
  if (!_track_map)
  {
    cerr << PHWHERE << " ERROR: Can't find " << _track_map_name.c_str() << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _svtx_seed_map = findNode::getClass<TrackSeedContainer>(topNode, "SvtxTrackSeedContainer");
  if (!_svtx_seed_map)
  {
    std::cout << "Creating node SvtxTrackSeedContainer" << std::endl;
    /// Get the DST Node
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

    /// Check that it is there
    if (!dstNode)
    {
      std::cerr << "DST Node missing, quitting" << std::endl;
      throw std::runtime_error("failed to find DST node in PHActsSourceLinks::createNodes");
    }

    /// Get the tracking subnode
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "SVTX"));

    /// Check that it is there
    if (!svtxNode)
    {
      svtxNode = new PHCompositeNode("SVTX");
      dstNode->addNode(svtxNode);
    }

    _svtx_seed_map = new TrackSeedContainer_v1();
    PHIODataNode<PHObject> *node = new PHIODataNode<PHObject>(_svtx_seed_map, "SvtxTrackSeedContainer", "PHObject");
    svtxNode->addNode(node);
  }

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, _cluster_map_name);
  if (!_cluster_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find node " <<_cluster_map_name << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!_tGeometry)
  {
    std::cout << PHWHERE << "Error, can't find acts tracking geometry" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

