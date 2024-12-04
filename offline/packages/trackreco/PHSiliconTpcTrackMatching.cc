#include "PHSiliconTpcTrackMatching.h"

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
PHSiliconTpcTrackMatching::PHSiliconTpcTrackMatching(const std::string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
{
  InitializeParameters();
}

//____________________________________________________________________________..
PHSiliconTpcTrackMatching::~PHSiliconTpcTrackMatching() = default;

//____________________________________________________________________________..
int PHSiliconTpcTrackMatching::InitRun(PHCompositeNode *topNode)
{
  UpdateParametersWithMacro();
  if(_test_windows)
  {
  _file = new TFile("track_match.root", "RECREATE");
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

  return ret;
}

//_____________________________________________________________________
void PHSiliconTpcTrackMatching::SetDefaultParameters()
{
  // Data on gasses @20 C and 760 Torr from the following source:
  // http://www.slac.stanford.edu/pubs/icfa/summer98/paper3/paper3.pdf
  // diffusion and drift velocity for 400kV for NeCF4 50/50 from calculations:
  // http://skipper.physics.sunysb.edu/~prakhar/tpc/HTML_Gases/split.html

  return;
}

//____________________________________________________________________________..
int PHSiliconTpcTrackMatching::process_event(PHCompositeNode * /*unused*/)
{
  if(Verbosity() > 2)
  {
    std::cout << " FIXME PHSiliconTpcTrackMatching " 
      << ( _zero_field ? "zero field is ON" : " zero field is OFF") << std::endl;
  }
  // _track_map contains the TPC seed track stubs
  // _track_map_silicon contains the silicon seed track stubs
  // _svtx_seed_map contains the combined silicon and tpc track seeds

  if (Verbosity() > 0)
  {
    cout << PHWHERE << " TPC track map size " << _track_map->size() << " Silicon track map size " << _track_map_silicon->size() << endl;
  }

  if (_track_map->size() == 0)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

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
  }

  // Find all matches of tpc and si tracklets in eta and phi, x and y
  //     If _pp_mode is not set, a match in z is also required - gives same behavior as old code
  std::multimap<unsigned int, unsigned int> tpc_matches;
  std::set<unsigned int> tpc_matched_set;
  std::set<unsigned int> tpc_unmatched_set;
  findEtaPhiMatches(tpc_matched_set, tpc_unmatched_set, tpc_matches);

  // Check that the crossing number is consistent with the tracklet z mismatch, removethe match otherwise
  // This does nothing if the crossing number is not set
  checkCrossingMatches(tpc_matches);

  // We have a complete list of all eta/phi matched tracks in the map "tpc_matches"
  // make the combined track seeds from tpc_matches
  for (auto [tpcid, si_id] : tpc_matches)
  {
    auto svtxseed = std::make_unique<SvtxTrackSeed_v2>();
    svtxseed->set_silicon_seed_index(si_id);
    svtxseed->set_tpc_seed_index(tpcid);
    // In pp mode, if a matched track does not have INTT clusters we have to find the crossing geometrically
    // Record the geometrically estimated crossing in the track seeds for later use if needed
    short int crossing_estimate = findCrossingGeometrically(tpcid, si_id);
    svtxseed->set_crossing_estimate(crossing_estimate);
    _svtx_seed_map->insert(svtxseed.get());

    if (Verbosity() > 1)
    {
      std::cout << "  combined seed id " << _svtx_seed_map->size() - 1 << " si id " << si_id << " tpc id " << tpcid << " crossing estimate " << crossing_estimate << std::endl;
    }
  }

  // Also make the unmatched TPC seeds into SvtxTrackSeeds
  for (auto tpcid : tpc_unmatched_set)
  {
    auto svtxseed = std::make_unique<SvtxTrackSeed_v2>();
    svtxseed->set_tpc_seed_index(tpcid);
    _svtx_seed_map->insert(svtxseed.get());

    if (Verbosity() > 1)
    {
      std::cout << "  converted unmatched TPC seed id " << _svtx_seed_map->size() - 1 << " tpc id " << tpcid << std::endl;
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
    }

    cout << "PHSiliconTpcTrackMatching::process_event(PHCompositeNode *topNode) Leaving process_event" << endl;
  }
  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}

short int PHSiliconTpcTrackMatching::findCrossingGeometrically(unsigned int tpcid, unsigned int si_id)
{
  // loop over all matches and check for ones with no INTT clusters in the silicon seed
  TrackSeed *si_track = _track_map_silicon->get(si_id);
  const short int crossing = si_track->get_crossing();
  const double si_z = TrackSeedHelper::get_z(si_track);

  TrackSeed *tpc_track = _track_map->get(tpcid);
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

double PHSiliconTpcTrackMatching::getBunchCrossing(unsigned int trid, double z_mismatch)
{
  double vdrift = _tGeometry->get_drift_velocity();  // cm/ns
  vdrift *= 1000.0;                                  // cm/microsecond
  //  double vdrift = 8.00;  // cm /microsecond
  // double z_bunch_separation = 0.106 * vdrift;  // 106 ns bunch crossing interval, as in pileup generator
  double z_bunch_separation = (crossing_period / 1000.0) * vdrift;  // 106 ns bunch crossing interval, as in pileup generator

  // The sign of z_mismatch will depend on which side of the TPC the tracklet is in
  TrackSeed *track = _track_map->get(trid);

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

int PHSiliconTpcTrackMatching::End(PHCompositeNode * /*unused*/)
{
  if(_test_windows)
  {
  _file->cd();
  _tree->Write();
  _file->Close();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSiliconTpcTrackMatching::GetNodes(PHCompositeNode *topNode)
{
  //---------------------------------
  // Get additional objects off the Node Tree
  //---------------------------------

  _cluster_crossing_map = findNode::getClass<TrkrClusterCrossingAssoc>(topNode, "TRKR_CLUSTERCROSSINGASSOC");
  if (!_cluster_crossing_map)
  {
    cerr << PHWHERE << " ERROR: Can't find TRKR_CLUSTERCROSSINGASSOC " << endl;
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

void PHSiliconTpcTrackMatching::findEtaPhiMatches(
    std::set<unsigned int> &tpc_matched_set,
    std::set<unsigned int> &tpc_unmatched_set,
    std::multimap<unsigned int, unsigned int> &tpc_matches)
{
  // loop over the TPC track seeds
  for (unsigned int phtrk_iter = 0;
       phtrk_iter < _track_map->size();
       ++phtrk_iter)
  {
    _tracklet_tpc = _track_map->get(phtrk_iter);
    if (!_tracklet_tpc)
    {
      continue;
    }

    unsigned int tpcid = phtrk_iter;
    if (Verbosity() > 1)
    {
      std::cout
          << __LINE__
          << ": Processing seed itrack: " << tpcid
          << ": nhits: " << _tracklet_tpc->size_cluster_keys()
          << ": Total tracks: " << _track_map->size()
          << ": phi: " << _tracklet_tpc->get_phi()
          << endl;
    }

    double tpc_phi, tpc_eta, tpc_pt;
    float tpc_px, tpc_py, tpc_pz;
    int tpc_q;
    Acts::Vector3 tpc_pos;
    if (_zero_field) {
      auto cluster_list = getTrackletClusterList(_tracklet_tpc);

      Acts::Vector3  mom;
      bool ok_track;

      std::tie(ok_track, tpc_phi, tpc_eta, tpc_pt, tpc_pos, mom) = 
        TrackFitUtils::zero_field_track_params(_tGeometry, _cluster_map, cluster_list);
      if (!ok_track) { continue; }
      tpc_px = mom.x();
      tpc_py = mom.y();
      tpc_pz = mom.z();
      tpc_q = -100;
    } else {
      tpc_phi = _tracklet_tpc->get_phi();
      tpc_eta = _tracklet_tpc->get_eta();
      tpc_pt = fabs(1. / _tracklet_tpc->get_qOverR()) * (0.3 / 100.) * fieldstrength;

      tpc_pos = TrackSeedHelper::get_xyz(_tracklet_tpc);

      tpc_px = _tracklet_tpc->get_px();
      tpc_py = _tracklet_tpc->get_py();
      tpc_pz = _tracklet_tpc->get_pz();

      tpc_q = _tracklet_tpc->get_charge();
    }
    double mag = getMatchingInflationFactor(tpc_pt);

    if (Verbosity() > 8)
    {
      std::cout << " tpc stub: " << tpcid << " eta " << tpc_eta << " phi " << tpc_phi << " pt " << tpc_pt << " tpc z " << TrackSeedHelper::get_z(_tracklet_tpc) << std::endl;
    }

    // this factor will increase the window size at low pT
    // otherwise the matching efficiency drops off at low pT


    if (_use_old_matching)  // for testing only
    {
      mag = 1.0;
      if (tpc_pt < 6.0)
      {
        mag = 2;
      }
      if (tpc_pt < 3.0)
      {
        mag = 4.0;
      }
      if (tpc_pt < 1.5)
      {
        mag = 6.0;
      }
      mag = 1.0;
    }

    if (Verbosity() > 3)
    {
      cout << "TPC tracklet:" << endl;
      _tracklet_tpc->identify();
    }


    bool matched = false;

    // Now search the silicon track list for a match in eta and phi
    for (unsigned int phtrk_iter_si = 0;
         phtrk_iter_si < _track_map_silicon->size();
         ++phtrk_iter_si)
    {
      _tracklet_si = _track_map_silicon->get(phtrk_iter_si);
      if (!_tracklet_si)
      {

        continue;
      }
      bool eta_match = false;

    double si_phi, si_eta, si_pt;
    float si_px, si_py, si_pz;
    int si_q;
    Acts::Vector3 si_pos;
    if (_zero_field) {
      auto cluster_list = getTrackletClusterList(_tracklet_si);

      Acts::Vector3  mom;
      bool ok_track;

      std::tie(ok_track, si_phi, si_eta, si_pt, si_pos, mom) = 
        TrackFitUtils::zero_field_track_params(_tGeometry, _cluster_map, cluster_list);
      if (!ok_track) { continue; }
      si_px = mom.x();
      si_py = mom.y();
      si_pz = mom.z();
      si_q = -100;
    } else {
      si_eta = _tracklet_si->get_eta();
      si_phi = _tracklet_si->get_phi();

      si_pos = TrackSeedHelper::get_xyz(_tracklet_si);
      si_px = _tracklet_si->get_px();
      si_py = _tracklet_si->get_py();
      si_pz = _tracklet_si->get_pz();
	    si_q = _tracklet_si->get_charge();
    }
	  int si_crossing = _tracklet_si->get_crossing();
    unsigned int siid = phtrk_iter_si;

  if(_test_windows)
  {
    float data[] = {
      (float) m_event, (float) si_crossing,
      (float) si_q, (float) si_phi, (float) si_eta, (float) si_pos.x(), (float) si_pos.y(), (float) si_pos.z(), (float) si_px, (float) si_py, (float) si_pz,
      (float) tpc_q, (float) tpc_phi, (float) tpc_eta, (float) tpc_pos.x(), (float) tpc_pos.y(), (float) tpc_pos.z(), (float) tpc_px, (float) tpc_py, (float) tpc_pz,
      (float) tpcid, (float) siid};
    _tree->Fill(data);
  }

      if (fabs(tpc_eta - si_eta) < _eta_search_win * mag)
      {
        eta_match = true;
      }
      if (!eta_match)
      {
        continue;
      }

      bool position_match = false;
      if (_pp_mode)
      {
        if( std::abs(tpc_pos.x() - si_pos.x()) < _x_search_win * mag && std::abs(tpc_pos.y() - si_pos.y()) < _y_search_win * mag)
        {
          position_match = true;
        }
      }
      else
      {
        if (
            fabs(tpc_pos.x() - si_pos.x()) < _x_search_win * mag && fabs(tpc_pos.y() - si_pos.y()) < _y_search_win * mag && fabs(tpc_pos.z() - si_pos.z()) < _z_search_win * mag)
        {
          position_match = true;
        }
      }

      if (!position_match)
      {
        continue;
      }

      bool phi_match = false;
      if (fabs(tpc_phi - si_phi) < _phi_search_win * mag)
      {
        phi_match = true;
      }
      if (fabs(fabs(tpc_phi - si_phi) - 2.0 * M_PI) < _phi_search_win * mag)
      {
        phi_match = true;
      }
      if (!phi_match)
      {
        continue;
      }
      if (Verbosity() > 3)
      {
        cout << " testing for a match for TPC track " << tpcid << " with pT " << _tracklet_tpc->get_pt()
             << " and eta " << _tracklet_tpc->get_eta() << " with Si track " << siid << " with crossing " << _tracklet_si->get_crossing() << endl;
        cout << " tpc_phi " << tpc_phi << " si_phi " << si_phi << " dphi " << tpc_phi - si_phi << " phi search " << _phi_search_win * mag << " tpc_eta " << tpc_eta
             << " si_eta " << si_eta << " deta " << tpc_eta - si_eta << " eta search " << _eta_search_win * mag << endl;
        std::cout << "      tpc x " << tpc_pos.x() << " si x " << si_pos.x() << " tpc y " << tpc_pos.y() << " si y " << si_pos.y() << " tpc_z " << tpc_pos.z() << " si z " << si_pos.z() << std::endl;
        std::cout << "      x search " << _x_search_win * mag << " y search " << _y_search_win * mag << " z search " << _z_search_win * mag << std::endl;
      }

      // got a match, add to the list
      // These stubs are matched in eta, phi, x and y already
      matched = true;
      tpc_matches.insert(std::make_pair(tpcid, siid));
      tpc_matched_set.insert(tpcid);

      if (Verbosity() > 1)
      {
        cout << " found a match for TPC track " << tpcid << " with Si track " << siid << endl;
        cout << "          tpc_phi " << tpc_phi << " si_phi " << si_phi << " phi_match " << phi_match
             << " tpc_eta " << tpc_eta << " si_eta " << si_eta << " eta_match " << eta_match << endl;
        std::cout << "      tpc x " << tpc_pos.x() << " si x " << si_pos.x() << " tpc y " << tpc_pos.y() << " si y " << si_pos.y() << " tpc_z " << tpc_pos.z() << " si z " << si_pos.z() << std::endl;
      }

      // temporary!
      if (_test_windows)
      {
        cout << " Try_silicon: crossing" << si_crossing <<  "  pt " << tpc_pt << " tpc_phi " << tpc_phi << " si_phi " << si_phi << " dphi " << tpc_phi - si_phi <<  "   si_q" << si_q << "   tpc_q" << tpc_q
             << " tpc_eta " << tpc_eta << " si_eta " << si_eta << " deta " << tpc_eta - si_eta << " tpc_x " << tpc_pos.x() << " tpc_y " << tpc_pos.y() << " tpc_z " << tpc_pos.z()
             << " dx " << tpc_pos.x() - si_pos.x() << " dy " << tpc_pos.y() - si_pos.y() << " dz " << tpc_pos.z() - si_pos.z()
			 << endl;
      }

    }
    // if no match found, keep tpc seed for fitting
    if (!matched)
    {
      if (Verbosity() > 1)
      {
        cout << "inserted unmatched tpc seed " << tpcid << endl;
      }
      tpc_unmatched_set.insert(tpcid);
    }
  }

  return;
}

void PHSiliconTpcTrackMatching::checkCrossingMatches(std::multimap<unsigned int, unsigned int> &tpc_matches)
{
  // if the  crossing was assigned correctly, the (crossing corrected) track position should satisfy the Z matching cut
  // this is a rough check that this is the case

  float vdrift = _tGeometry->get_drift_velocity();

  std::multimap<unsigned int, unsigned int> bad_map;

  for (auto [tpcid, si_id] : tpc_matches)
  {
    TrackSeed *tpc_track = _track_map->get(tpcid);
    TrackSeed *si_track = _track_map_silicon->get(si_id);
    short int crossing = si_track->get_crossing();

    if (crossing == SHRT_MAX)
    {
      if (Verbosity() > 2)
      {
        std::cout << " drop si_track " << si_id << " with eta " << si_track->get_eta() << " and z " << TrackSeedHelper::get_z(si_track) << " because crossing is undefined " << std::endl;
      }
      continue;
    }

    float z_si = TrackSeedHelper::get_z(si_track);
    float z_tpc = TrackSeedHelper::get_z(tpc_track);
    float z_mismatch = z_tpc - z_si;

    float mag_crossing_z_mismatch = fabs(crossing) * crossing_period * vdrift;

    // We do not know the sign  of the z mismatch for a given crossing unless we know the drift direction in the TPC, use magnitude
    // could instead look up any TPC cluster key in the track to get side
    // z-mismatch can occasionally be up to 2 crossings due to TPC extrapolation precision
    if (fabs(fabs(z_mismatch) - mag_crossing_z_mismatch) < 3.0)
    {
      if (Verbosity() > 1)
      {
        std::cout << "  Success:  crossing " << crossing << " tpcid " << tpcid << " si id " << si_id
                  << " tpc z " << z_tpc << " si z " << z_si << " z_mismatch " << z_mismatch
                  << " mag_crossing_z_mismatch " << mag_crossing_z_mismatch << " drift velocity " << vdrift << std::endl;
      }
    }
    else
    {
      if (Verbosity() > 1)
      {
        std::cout << "  FAILURE:  crossing " << crossing << " tpcid " << tpcid << " si id " << si_id
                  << " tpc z " << z_tpc << " si z " << z_si << " z_mismatch " << z_mismatch
                  << " mag_crossing_z_mismatch " << mag_crossing_z_mismatch << std::endl;
      }

      // bad_map.insert(std::make_pair(tpcid, si_id));
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

double PHSiliconTpcTrackMatching::getMatchingInflationFactor(double tpc_pt)
{
  double mag = 1.0;

  if (tpc_pt > _match_function_ptmin)
  {
    mag = _match_function_a + _match_function_b / pow(tpc_pt, _match_function_pow);
  }

  // std::cout << "  tpc_pt = " << tpc_pt << " mag " << mag << " a " << _match_function_a << " b " << _match_function_b << std::endl;

  return mag;
}

std::vector<TrkrDefs::cluskey> PHSiliconTpcTrackMatching::getTrackletClusterList(TrackSeed* tracklet)
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

    // drop INTT clusters for now  -- TEMPORARY!
    // Note: the zerofield fit uses the INTT for xy fit but not yz fit
    /* if (layer > 2 && layer < 7) */
    /* { */
      /* continue; */
    /* } */


    cluskey_vec.push_back(key);
  }  // end loop over clusters for this track
  return cluskey_vec;
}
