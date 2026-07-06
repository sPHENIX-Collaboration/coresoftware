#include "PHSiliconTpcTrackMatchingDummy.h"

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
PHSiliconTpcTrackMatchingDummy::PHSiliconTpcTrackMatchingDummy(const std::string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
{
  InitializeParameters();
}

//____________________________________________________________________________..
PHSiliconTpcTrackMatchingDummy::~PHSiliconTpcTrackMatchingDummy() = default;

//____________________________________________________________________________..
int PHSiliconTpcTrackMatchingDummy::InitRun(PHCompositeNode *topNode)
{
  UpdateParametersWithMacro();
  if(_test_windows)
  {
  _file = new TFile(_file_name.c_str(), "RECREATE");
  _tree = new TNtuple("track_match", "track_match",
                      "event:sicrossing:siq:siphi:sieta:six:siy:siz:sipx:sipy:sipz:tpcq:tpcphi:tpceta:tpcx:tpcy:tpcz:tpcpx:tpcpy:tpcpz:tpcid:siid");
  }
  // put these in the output file
  cout << PHWHERE << " Search windows: phi " << _phi_search_win << " eta "
       << _eta_search_win << " _pp_mode " << _pp_mode << " _use_intt_crossing " << _use_intt_crossing << endl;

  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }
  std::istringstream stringline(m_fieldMap);
  stringline >> fieldstrength;

  // initialize the WindowMatchers
  window_dx.init_bools("dx", _print_windows || Verbosity()   >0);
  window_dy.init_bools("dy", _print_windows || Verbosity()   >0);
  window_dz.init_bools("dz", _print_windows || Verbosity()   >0);
  window_dphi.init_bools("dphi", _print_windows || Verbosity() >0);
  window_deta.init_bools("deta", _print_windows || Verbosity() >0);


  return ret;
}

//_____________________________________________________________________
void PHSiliconTpcTrackMatchingDummy::SetDefaultParameters()
{
  // Data on gasses @20 C and 760 Torr from the following source:
  // http://www.slac.stanford.edu/pubs/icfa/summer98/paper3/paper3.pdf
  // diffusion and drift velocity for 400kV for NeCF4 50/50 from calculations:
  // http://skipper.physics.sunysb.edu/~prakhar/tpc/HTML_Gases/split.html

  return;
}

std::string PHSiliconTpcTrackMatchingDummy::WindowMatcher::print_fn(const Arr3D& dat) {
  std::ostringstream os;
  if (dat[1]==0.) {
    os << dat[0];
  } else {
    os << dat[0] << (dat[1]>0 ? "+" : "") << dat[1] <<"*exp("<<dat[2]<<"/pT)";
  }
  return os.str();
}

void PHSiliconTpcTrackMatchingDummy::WindowMatcher::init_bools(const std::string& tag, const bool print) {
  // set values for positive tracks
  fabs_max_posQ = (posLo[0]==100.);
  posLo_b0 = (posLo[1]==0.);
  posHi_b0 = (posHi[1]==0.);
  // if no values for negative tracks, copy over from positive tracks
  if (negHi[0]==100.) {
    negLo = posLo;
    negHi = posHi;
    fabs_max_negQ = fabs_max_posQ;
    negLo_b0 = posLo_b0;
    negHi_b0 = posHi_b0;
    min_pt_negQ = min_pt_posQ;
  } else {
    fabs_max_negQ = (negLo[0]==100.);
    negLo_b0 = (negLo[1]==0.);
    negHi_b0 = (negHi[1]==0.);
  }
  if (print) {
    std::cout << " Track matching window, " << tag << ":" << std::endl;

    if (posHi==negHi && posLo == negLo) {
      std::cout << "  all tracks: ";
    } else {
      std::cout << "   +Q tracks: ";
    }
    if (posLo[0]==100) {
      std::cout << "  |" << tag <<"| < " << print_fn(posHi) << std::endl;
    } else {
      std::cout << print_fn(posLo) <<" < " << tag << " < " << print_fn(posHi) << std::endl;
    }

    if (posHi != negHi || posLo != negLo) {
      std::cout << "   -Q tracks: ";
      if (negLo[0]==100) {
        std::cout << "  |" << tag <<"| < " << print_fn(negHi) << std::endl;
      } else {
        std::cout << print_fn(negLo) <<" < " << tag << " < " << print_fn(negHi) << std::endl;
      }
    }
  }

}

bool PHSiliconTpcTrackMatchingDummy::WindowMatcher::in_window
(const bool posQ, const double tpc_pt, const double tpc_X, const double si_X)
{
  const auto delta = tpc_X-si_X;
  
  if (posQ) {
    double pt = (tpc_pt<min_pt_posQ) ? min_pt_posQ : tpc_pt;
    if (fabs_max_posQ) {
      return fabs(delta) < fn_exp(posHi, posHi_b0, pt);
    } else {
      return (delta > fn_exp(posLo, posLo_b0, pt)
           && delta < fn_exp(posHi, posHi_b0, pt));
    }
  } else {
    double pt = (tpc_pt<min_pt_negQ) ? min_pt_negQ : tpc_pt;
    if (fabs_max_negQ) {
      return fabs(delta) < fn_exp(negHi, negHi_b0, pt);
    } else {
      return (delta > fn_exp(negLo, negLo_b0, pt)
           && delta < fn_exp(negHi, negHi_b0, pt));
    }
  }
}

//____________________________________________________________________________..
int PHSiliconTpcTrackMatchingDummy::process_event(PHCompositeNode * /*unused*/)
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

  if (_track_map->size() == 0)
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

    
    auto svtxseed = std::make_unique<SvtxTrackSeed_v2>();
    svtxseed->set_silicon_seed_index(si_id);
    svtxseed->set_tpc_seed_index(si_id);
    // In pp mode, if a matched track does not have INTT clusters we have to find the crossing geometrically
    // Record the geometrically estimated crossing in the track seeds for later use if needed
    svtxseed->set_crossing_estimate(crossing);
    _track_map->insert(dummy.get());
    _svtx_seed_map->insert(svtxseed.get());

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

    cout << "PHSiliconTpcTrackMatchingDummy::process_event(PHCompositeNode *topNode) Leaving process_event" << endl;
  }
  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}

double PHSiliconTpcTrackMatchingDummy::getBunchCrossing(unsigned int trid, double z_mismatch)
{
  const double vdrift = _tGeometry->get_drift_velocity();  // cm/ns
  const double z_bunch_separation = sphenix_constants::time_between_crossings * vdrift; // cm

  // The sign of z_mismatch will depend on which side of the TPC the tracklet is in
  TrackSeed *track = _track_map->get(trid);

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

int PHSiliconTpcTrackMatchingDummy::End(PHCompositeNode * /*unused*/)
{
  if(_test_windows)
  {
  _file->cd();
  _tree->Write();
  _file->Close();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSiliconTpcTrackMatchingDummy::GetNodes(PHCompositeNode *topNode)
{
  //---------------------------------
  // Get additional objects off the Node Tree
  //---------------------------------

  _cluster_crossing_map = findNode::getClass<TrkrClusterCrossingAssoc>(topNode, "TRKR_CLUSTERCROSSINGASSOC");
  if (!_cluster_crossing_map)
  {
    //cerr << PHWHERE << " ERROR: Can't find TRKR_CLUSTERCROSSINGASSOC " << endl;
    // return Fun4AllReturnCodes::ABORTEVENT;
  }

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

void PHSiliconTpcTrackMatchingDummy::checkZMatches(
    std::multimap<unsigned int, unsigned int> &tpc_matches,
    std::multimap<unsigned int, unsigned int> &bad_map)
{
  // for _pp_mode=false, assume zero crossings for all track matches
  // for _pp_mode=true, do crossing correction on track position z according to side and vdrift
  // z matching criteria follows window_z
  // there is a dz threshold cut to avoid window_z blow up at low pT

  float vdrift = _tGeometry->get_drift_velocity();

  for (auto [tpcid, si_id] : tpc_matches)
  {
    TrackSeed *tpc_track = _track_map->get(tpcid);
    TrackSeed *si_track = _track_map_silicon->get(si_id);

    short int crossing = si_track->get_crossing();
    float tpc_pt, tpc_z, si_z;
    int tpc_q;
    if (_zero_field) {
      auto cluster_list_tpc = getTrackletClusterList(tpc_track);
      auto cluster_list_si = getTrackletClusterList(si_track);

      tpc_pt = std::get<3>(TrackFitUtils::zero_field_track_params(_tGeometry, _cluster_map, cluster_list_tpc));
      tpc_z = std::get<4>(TrackFitUtils::zero_field_track_params(_tGeometry, _cluster_map, cluster_list_tpc)).z();
      tpc_q = -100;

      si_z = std::get<4>(TrackFitUtils::zero_field_track_params(_tGeometry, _cluster_map, cluster_list_si)).z();
    } else {
      tpc_pt = fabs(1. / _tracklet_tpc->get_qOverR()) * (0.3 / 100.) * fieldstrength;
      tpc_z = TrackSeedHelper::get_z(tpc_track);
      tpc_q = _tracklet_tpc->get_charge();
      si_z = TrackSeedHelper::get_z(si_track);
    }

    // get TPC side from one of the TPC clusters
    std::vector<TrkrDefs::cluskey> temp_clusters = getTrackletClusterList(tpc_track);
    if(temp_clusters.size() == 0) { continue; }
    unsigned int this_side = TpcDefs::getSide(temp_clusters[0]);

    bool is_posQ = (tpc_q>0.);

    float z_mismatch = tpc_z - si_z;
    float tpc_z_corrected = _clusterCrossingCorrection.correctZ(tpc_z, this_side, crossing);
    float z_mismatch_corrected = tpc_z_corrected - si_z;

    bool z_match = false;
    if (_pp_mode)
    {
      if (crossing == SHRT_MAX)
      {
        if (Verbosity() > 2)
        {
          std::cout << " drop si_track " << si_id << " with eta " << si_track->get_eta() << " and z " << TrackSeedHelper::get_z(si_track) << " because crossing is undefined " << std::endl;
        }
        continue;
      }

      if (window_dz.in_window(is_posQ, tpc_pt, tpc_z_corrected, si_z) && (fabs(z_mismatch_corrected) < _crossing_deltaz_max))
      {
        z_match = true;
      }
      else if (fabs(z_mismatch_corrected) < _crossing_deltaz_min)
      {
	z_match = true;
      }
    }
    else
    {
      if (window_dz.in_window(is_posQ, tpc_pt, tpc_z, si_z) && (fabs(z_mismatch) < _crossing_deltaz_max))
      {
        z_match = true;
      }
      else if (fabs(z_mismatch) < _crossing_deltaz_min)
      {
	z_match = true;
      }
    }

    if (z_match)
    {
      if (Verbosity() > 1)
      {
        std::cout << "  Success:  crossing " << crossing << " tpcid " << tpcid << " si id " << si_id
                  << " tpc z " << tpc_z << " si z " << si_z << " z_mismatch " << z_mismatch << "tpc z corrected " << tpc_z_corrected
                  << " z_mismatch_corrected " << z_mismatch_corrected << " drift velocity " << vdrift << std::endl;
      }
    }
    else
    {
      if (Verbosity() > 1)
      {
        std::cout << "  FAILURE:  crossing " << crossing << " tpcid " << tpcid << " si id " << si_id
                  << " tpc z " << tpc_z << " si z " << si_z << " z_mismatch " << z_mismatch << "tpc_z_corrected " << tpc_z_corrected
                  << " z_mismatch_corrected " << z_mismatch_corrected << std::endl;
      }

      bad_map.insert(std::make_pair(tpcid, si_id));
    }
  }

  // remove bad entries from tpc_matches
  for (auto [tpcid, si_id] : bad_map)
  {
    // Have to iterate over tpc_matches and examine each pair to find the one matching bad_map
    // this logic works because we call the equal range on vertex_map for every id_pair
    // so we only delete one entry per equal range call
    auto ret = tpc_matches.equal_range(tpcid);
    for (auto it = ret.first; it != ret.second; ++it)
    {
      if (it->first == tpcid && it->second == si_id)
      {
        if (Verbosity() > 1)
        {
          std::cout << "                        erasing tpc_matches entry for tpcid " << tpcid << " si_id " << si_id << std::endl;
        }
        tpc_matches.erase(it);
        break;  // the iterator is no longer valid
      }
    }
  }

  return;
}

std::vector<TrkrDefs::cluskey> PHSiliconTpcTrackMatchingDummy::getTrackletClusterList(TrackSeed* tracklet)
{
  std::vector<TrkrDefs::cluskey> cluskey_vec;
  for (auto clusIter = tracklet->begin_cluster_keys();
       clusIter != tracklet->end_cluster_keys();
       ++clusIter)
  {
    auto key = *clusIter;
    auto cluster = _cluster_map->findCluster(key);
    if (!cluster)
    {
      if(Verbosity() > 1)
      {
        std::cout << PHWHERE << "Failed to get cluster with key " << key << std::endl;
      }
      continue;
    }

    /// Make a safety check for clusters that couldn't be attached to a surface
    auto surf = _tGeometry->maps().getSurface(key, cluster);
    if (!surf)
    {
      continue;
    }

    // drop some bad layers in the TPC completely
    unsigned int layer = TrkrDefs::getLayer(key);
    if (layer == 7 || layer == 22 || layer == 23 || layer == 38 || layer == 39)
    {
      continue;
    }

    cluskey_vec.push_back(key);
  }  // end loop over clusters for this track
  return cluskey_vec;
}
