#include "PHTrackPruner.h"

/// Tracking includes
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrDefs.h>  // for cluskey, getTrkrId, tpcId

#include <trackbase_historic/SvtxTrackSeed_v2.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v2.h>
#include <trackbase_historic/TrackSeedHelper.h>

#include <trackbase_historic/SvtxTrack.h>  // for SvtxTrack
#include <trackbase_historic/SvtxTrackMap.h>

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

namespace
{
  //! get cluster keys from a given track
  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track)
  {
    std::vector<TrkrDefs::cluskey> out;
    for (const auto& seed : {track->get_silicon_seed(), track->get_tpc_seed()})
    {
      if (seed)
      {
        std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
      }
    }
    return out;
  }

  //! get cluster keys from a given track associated with states
  std::vector<TrkrDefs::cluskey> get_state_keys(SvtxTrack* track)
  {
    std::vector<TrkrDefs::cluskey> out;
    for (auto state_iter = track->begin_states();
         state_iter != track->end_states();
         ++state_iter)
    {
      SvtxTrackState* tstate = state_iter->second;
      auto stateckey = tstate->get_cluskey();
      out.push_back(stateckey);
    }
    return out;
  }

  /// return number of clusters of a given type that belong to a tracks
  template <int type>
  int count_clusters(const std::vector<TrkrDefs::cluskey>& keys)
  {
    return std::count_if(keys.begin(), keys.end(),
                         [](const TrkrDefs::cluskey& key)
                         { return TrkrDefs::getTrkrId(key) == type; });
  }
}

//____________________________________________________________________________..
PHTrackPruner::PHTrackPruner(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
PHTrackPruner::~PHTrackPruner() = default;

//____________________________________________________________________________..
int PHTrackPruner::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  return ret;
}

//____________________________________________________________________________..
int PHTrackPruner::process_event(PHCompositeNode * /*unused*/)
{
  // _tpc_seed_map contains the TPC seed track stubs
  // _si_seed_map contains the silicon seed track stubs
  // _svtx_seed_map contains the combined silicon and tpc track seeds
  // _svtx_track_map contains the fitted acts track stubs

  if (Verbosity() > 0)
  {
    cout << PHWHERE << " TPC seed map size " << _tpc_seed_map->size()
	 << " Silicon seed map size " << _si_seed_map->size()
	 << " Svtx seed map size " << _svtx_seed_map->size()
	 << " Svtx track map size " << _svtx_track_map->size()
	 << endl;
  }

  if (_svtx_track_map->size() == 0)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  std::multimap<unsigned int, unsigned int> good_matches;

  for (auto &iter : *_svtx_track_map)
  {
    _svtx_track = iter.second;

    if(!checkTrack(_svtx_track))
    {
      continue;
    }
    if (Verbosity() > 1) { std::cout<<"Pass track selection"<<std::endl; }

    _tpc_seed = _svtx_track->get_tpc_seed();
    _si_seed = _svtx_track->get_silicon_seed();
    if (_tpc_seed && _si_seed)
    {
      if (Verbosity() > 1) { std::cout<<"Insert tpcid and siid into good_matches"<<std::endl; }
      int tpcid = _tpc_seed_map->find(_tpc_seed);
      int siid = _si_seed_map->find(_si_seed);
      good_matches.insert(std::make_pair(tpcid, siid));
    }
  }

  for (auto [tpcid, siid] : good_matches)
  {
      if (Verbosity() > 1) { std::cout<<"Insert pruned svtx seed map"<<std::endl; }
    auto _svtx_seed = std::make_unique<SvtxTrackSeed_v2>();
    _svtx_seed->set_silicon_seed_index(siid);
    _svtx_seed->set_tpc_seed_index(tpcid);
    // In pp mode, if a matched track does not have INTT clusters we have to find the crossing geometrically
    // Record the geometrically estimated crossing in the track seeds for later use if needed
    short int crossing_estimate = findCrossingGeometrically(tpcid, siid);
    _svtx_seed->set_crossing_estimate(crossing_estimate);
    _pruned_svtx_seed_map->insert(_svtx_seed.get());

    if (Verbosity() > 1)
    {
      std::cout << "  combined seed id " << _pruned_svtx_seed_map->size() - 1 << " si id " << siid << " tpc id " << tpcid << " crossing estimate " << crossing_estimate << std::endl;
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << "final svtx seed map size " << _pruned_svtx_seed_map->size() << std::endl;
  }

  if (Verbosity() > 1)
  {
    for (const auto &seed : *_pruned_svtx_seed_map)
    {
      seed->identify();
    }

    cout << "PHTrackPruner::process_event(PHCompositeNode *topNode) Leaving process_event" << endl;
  }
  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}

bool PHTrackPruner::checkTrack(SvtxTrack *track)
{
  if(!track)
  {
    if (Verbosity() > 1) { std::cout<<"invalid track"<<std::endl; }
    return false;
  }

  //pt cut
  if(track->get_pt() < m_track_pt_low_cut)
  {
    if (Verbosity() > 1) { std::cout<<"Track pt "<<track->get_pt()<<" , pt cut "<<m_track_pt_low_cut<<std::endl; }
    return false;
  }

  //quality cut
  if(track->get_quality() > m_track_quality_high_cut)
  {
    if (Verbosity() > 1) { std::cout<<"Track quality "<<track->get_quality()<<" , quality cut "<<m_track_quality_high_cut<<std::endl; }
    return false;
  }

  //number of clusters cut
  const auto cluster_keys(get_cluster_keys(track));
  if (count_clusters<TrkrDefs::mvtxId>(cluster_keys) < m_nmvtx_clus_low_cut)
  {
    if (Verbosity() > 1) { std::cout<<"nmvtx "<<count_clusters<TrkrDefs::mvtxId>(cluster_keys)<<" , nmvtx cut "<<m_nmvtx_clus_low_cut<<std::endl; }
    return false;
  }
  if (count_clusters<TrkrDefs::inttId>(cluster_keys) < m_nintt_clus_low_cut)
  {
    if (Verbosity() > 1) { std::cout<<"nintt "<<count_clusters<TrkrDefs::inttId>(cluster_keys)<<" , nintt cut "<<m_nintt_clus_low_cut<<std::endl; }
    return false;
  }
  if (count_clusters<TrkrDefs::tpcId>(cluster_keys) < m_ntpc_clus_low_cut)
  {
    if (Verbosity() > 1) { std::cout<<"ntpc "<<count_clusters<TrkrDefs::tpcId>(cluster_keys)<<" , ntpc cut "<<m_ntpc_clus_low_cut<<std::endl; }
    return false;
  }
  if (count_clusters<TrkrDefs::micromegasId>(cluster_keys) < m_ntpot_clus_low_cut)
  {
    if (Verbosity() > 1) { std::cout<<"nmicromegas "<<count_clusters<TrkrDefs::micromegasId>(cluster_keys)<<" , nmicromegas cut "<<m_ntpot_clus_low_cut<<std::endl; }
    return false;
  }

  //number of states cut
  const auto state_keys(get_state_keys(track));
  if (count_clusters<TrkrDefs::mvtxId>(state_keys) < m_nmvtx_states_low_cut)
  {
    if (Verbosity() > 1) { std::cout<<"nmvtxstates "<<count_clusters<TrkrDefs::mvtxId>(state_keys)<<" , nmvtxstates cut "<<m_nmvtx_states_low_cut<<std::endl; }
    return false;
  }
  if (count_clusters<TrkrDefs::inttId>(state_keys) < m_nintt_states_low_cut)
  {
    if (Verbosity() > 1) { std::cout<<"ninttstates "<<count_clusters<TrkrDefs::inttId>(state_keys)<<" , ninttstates cut "<<m_nintt_states_low_cut<<std::endl; }
    return false;
  }
  if (count_clusters<TrkrDefs::tpcId>(state_keys) < m_ntpc_states_low_cut)
  {
    if (Verbosity() > 1) { std::cout<<"ntpcstates "<<count_clusters<TrkrDefs::tpcId>(state_keys)<<" , ntpcstates cut "<<m_ntpc_states_low_cut<<std::endl; }
    return false;
  }
  if (count_clusters<TrkrDefs::micromegasId>(state_keys) < m_ntpot_states_low_cut)
  {
    if (Verbosity() > 1) { std::cout<<"nmicromegasstates "<<count_clusters<TrkrDefs::micromegasId>(state_keys)<<" , nmicromegasstates cut "<<m_ntpot_states_low_cut<<std::endl; }
    return false;
  }

  return true;
}

int PHTrackPruner::End(PHCompositeNode * /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackPruner::GetNodes(PHCompositeNode *topNode)
{
  //---------------------------------
  // Get additional objects off the Node Tree
  //---------------------------------

  _svtx_track_map = findNode::getClass<SvtxTrackMap>(topNode, _svtx_track_map_name);
  if (!_svtx_track_map)
  {
    cerr << PHWHERE << " ERROR: Can't find " << _svtx_track_map_name.c_str()  << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _si_seed_map = findNode::getClass<TrackSeedContainer>(topNode, _si_seed_map_name);
  if (!_si_seed_map)
  {
    cerr << PHWHERE << " ERROR: Can't find " << _si_seed_map_name.c_str()  << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _tpc_seed_map = findNode::getClass<TrackSeedContainer>(topNode, _tpc_seed_map_name);
  if (!_tpc_seed_map)
  {
    cerr << PHWHERE << " ERROR: Can't find " << _tpc_seed_map_name.c_str() << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _svtx_seed_map = findNode::getClass<TrackSeedContainer>(topNode, _svtx_seed_map_name);
  if (!_svtx_seed_map)
  {
    cerr << PHWHERE << " ERROR: Can't find " << _svtx_seed_map_name.c_str() << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _pruned_svtx_seed_map = findNode::getClass<TrackSeedContainer>(topNode, _pruned_svtx_seed_map_name);
  if (!_pruned_svtx_seed_map)
  {
    std::cout << "Creating node " << _pruned_svtx_seed_map_name.c_str() << std::endl;
    /// Get the DST Node
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

    /// Check that it is there
    if (!dstNode)
    {
      std::cerr << "DST Node missing, quitting" << std::endl;
      throw std::runtime_error("failed to find DST node in PHTrackPruner::GetNodes");
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

    _pruned_svtx_seed_map = new TrackSeedContainer_v1();
    PHIODataNode<PHObject> *node = new PHIODataNode<PHObject>(_pruned_svtx_seed_map, _pruned_svtx_seed_map_name, "PHObject");
    svtxNode->addNode(node);
  }

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
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

short int PHTrackPruner::findCrossingGeometrically(unsigned int tpcid, unsigned int si_id)
{
  // loop over all matches and check for ones with no INTT clusters in the silicon seed
  TrackSeed *si_track = _si_seed_map->get(si_id);
  const short int crossing = si_track->get_crossing();
  const double si_z = TrackSeedHelper::get_z(si_track);

  TrackSeed *tpc_track = _tpc_seed_map->get(tpcid);
  const double tpc_z = TrackSeedHelper::get_z(tpc_track);

  // this is an initial estimate of the bunch crossing based on the z-mismatch of the seeds for this track
  short int crossing_estimate = (short int) getBunchCrossing(tpcid, tpc_z - si_z);

  if (Verbosity() > 1)
  {
    std::cout << "findCrossing: "
              << " tpcid " << tpcid << " si_id " << si_id << " tpc_z " << tpc_z << " si_z " << si_z << " dz " << tpc_z - si_z
              << " INTT crossing " << crossing << " crossing_estimate " << crossing_estimate << std::endl;
  }

  return crossing_estimate;
}

double PHTrackPruner::getBunchCrossing(unsigned int trid, double z_mismatch)
{
  const double vdrift = _tGeometry->get_drift_velocity();  // cm/ns
  const double z_bunch_separation = sphenix_constants::time_between_crossings * vdrift; // cm

  // The sign of z_mismatch will depend on which side of the TPC the tracklet is in
  TrackSeed *track = _tpc_seed_map->get(trid);

  // crossing
  double crossings = z_mismatch / z_bunch_separation;

  // Check the TPC side for the first cluster in the track
  unsigned int side = 10;
  std::set<short int> side_set;
  for (TrackSeed::ConstClusterKeyIter iter = track->begin_cluster_keys();
       iter != track->end_cluster_keys();
       ++iter)
  {
    TrkrDefs::cluskey cluster_key = *iter;
    unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);
    if (trkrid == TrkrDefs::tpcId)
    {
      side = TpcDefs::getSide(cluster_key);
      side_set.insert(side);
    }
  }

  if (side == 10)
  {
    return SHRT_MAX;
  }

  if (side_set.size() == 2 && Verbosity() > 1)
  {
    std::cout << "     WARNING: tpc seed " << trid << " changed TPC sides, "
              << "  final side " << side << std::endl;
  }

  // if side = 1 (north, +ve z side), a positive t0 will make the cluster late relative to true z, so it will look like z is less positive
  // so a negative z mismatch for side 1 means a positive t0, and positive crossing, so reverse the sign for side 1
  if (side == 1)
  {
    crossings *= -1.0;
  }

  if (Verbosity() > 1)
  {
    std::cout << "  gettrackid " << trid << " side " << side << " z_mismatch " << z_mismatch << " crossings " << crossings << std::endl;
  }

  return crossings;
}
