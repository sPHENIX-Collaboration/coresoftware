#include "PHSiliconTpcTrackMatching.h"

/// Tracking includes
#include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterCrossingAssoc.h>
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrDefs.h>  // for cluskey, getTrkrId, tpcId

#include <trackbase_historic/SvtxTrackSeed_v2.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v2.h>

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

    // returns SHRT_MAX if no INTT clusters in silicon seed
    short int crossing = getCrossingIntt(_tracklet_si);
    if (_use_intt_crossing)
    {
      _tracklet_si->set_crossing(crossing);
    }  // flag is for testing only, use_intt_crossing should always be true!

    if (Verbosity() > 8)
    {
      std::cout << " silicon stub: " << trackid << " eta " << _tracklet_si->get_eta() << " pt " << _tracklet_si->get_pt() << " si z " << _tracklet_si->get_z() << " crossing " << crossing << std::endl;
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

  return Fun4AllReturnCodes::EVENT_OK;
}

short int PHSiliconTpcTrackMatching::findCrossingGeometrically(unsigned int tpcid, unsigned int si_id)
{
  // loop over all matches and check for ones with no INTT clusters in the silicon seed
  TrackSeed *si_track = _track_map_silicon->get(si_id);
  short int crossing = si_track->get_crossing();

  double si_z = si_track->get_z();
  TrackSeed *tpc_track = _track_map->get(tpcid);
  double tpc_z = tpc_track->get_z();

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

    double tpc_phi = _tracklet_tpc->get_phi();
    double tpc_eta = _tracklet_tpc->get_eta();
    double tpc_pt = fabs(1. / _tracklet_tpc->get_qOverR()) * (0.3 / 100.) * fieldstrength;
    if (Verbosity() > 8)
    {
      std::cout << " tpc stub: " << tpcid << " eta " << tpc_eta << " phi " << tpc_phi << " pt " << tpc_pt << " tpc z " << _tracklet_tpc->get_z() << std::endl;
    }

    // this factor will increase the window size at low pT
    // otherwise the matching efficiency drops off at low pT

    double mag = getMatchingInflationFactor(tpc_pt);

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
    }

    if (Verbosity() > 3)
    {
      cout << "TPC tracklet:" << endl;
      _tracklet_tpc->identify();
    }

    double tpc_x = _tracklet_tpc->get_x();
    double tpc_y = _tracklet_tpc->get_y();
    double tpc_z = _tracklet_tpc->get_z();

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
      double si_eta = _tracklet_si->get_eta();
      if (fabs(tpc_eta - si_eta) < _eta_search_win * mag)
      {
        eta_match = true;
      }
      if (!eta_match)
      {
        continue;
      }
      unsigned int siid = phtrk_iter_si;
      double si_x = _tracklet_si->get_x();
      double si_y = _tracklet_si->get_y();
      double si_z = _tracklet_si->get_z();
      bool position_match = false;
      if (_pp_mode)
      {
        if (
            fabs(tpc_x - si_x) < _x_search_win * mag && fabs(tpc_y - si_y) < _y_search_win * mag)
        {
          position_match = true;
        }
      }
      else
      {
        if (
            fabs(tpc_x - si_x) < _x_search_win * mag && fabs(tpc_y - si_y) < _y_search_win * mag && fabs(tpc_z - si_z) < _z_search_win * mag)
        {
          position_match = true;
        }
      }

      if (!position_match)
      {
        continue;
      }

      bool phi_match = false;
      double si_phi = _tracklet_si->get_phi();
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
        std::cout << "      tpc x " << tpc_x << " si x " << si_x << " tpc y " << tpc_y << " si y " << si_y << " tpc_z " << tpc_z << " si z " << si_z << std::endl;
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
        std::cout << "      tpc x " << tpc_x << " si x " << si_x << " tpc y " << tpc_y << " si y " << si_y << " tpc_z " << tpc_z << " si z " << si_z << std::endl;
      }

      // temporary!
      if (_test_windows)
      {
        cout << " Try_silicon:  pt " << tpc_pt << " tpc_phi " << tpc_phi << " si_phi " << si_phi << " dphi " << tpc_phi - si_phi
             << " tpc_eta " << tpc_eta << " si_eta " << si_eta << " deta " << tpc_eta - si_eta << " tpc_x " << tpc_x << " tpc_y " << tpc_y << " tpc_z " << tpc_z
             << " dx " << tpc_x - si_x << " dy " << tpc_y - si_y << " dz " << tpc_z - si_z
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

short int PHSiliconTpcTrackMatching::getCrossingIntt(TrackSeed *si_track)
{
  // If the Si track contains an INTT hit, use it to get the bunch crossing offset

  std::vector<short int> intt_crossings = getInttCrossings(si_track);

  bool keep_it = true;
  short int crossing_keep = 0;
  if (intt_crossings.size() == 0)
  {
    keep_it = false;
  }
  else
  {
    crossing_keep = intt_crossings[0];
    for (unsigned int ic = 1; ic < intt_crossings.size(); ++ic)
    {
      if (intt_crossings[ic] != crossing_keep)
      {
        if (Verbosity() > 1)
        {
          std::cout << " Warning: INTT crossings not all the same "
                    << " crossing_keep " << crossing_keep << " new crossing " << intt_crossings[ic] << " keep the first one in the list" << std::endl;
        }
      }
    }
  }

  if (keep_it)
  {
    return crossing_keep;
  }

  return SHRT_MAX;
}

std::vector<short int> PHSiliconTpcTrackMatching::getInttCrossings(TrackSeed *si_track)
{
  std::vector<short int> intt_crossings;

  // If the Si track contains an INTT hit, use it to get the bunch crossing offset
  // loop over associated clusters to get keys for silicon cluster
  for (TrackSeed::ConstClusterKeyIter iter = si_track->begin_cluster_keys();
       iter != si_track->end_cluster_keys();
       ++iter)
  {
    TrkrDefs::cluskey cluster_key = *iter;
    const unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);

    if (Verbosity() > 1)
    {
      unsigned int layer = TrkrDefs::getLayer(cluster_key);

      if (trkrid == TrkrDefs::mvtxId)
      {
        TrkrCluster *cluster = _cluster_map->findCluster(cluster_key);
        if (!cluster)
        {
          continue;
        }

        Acts::Vector3 global = _tGeometry->getGlobalPosition(cluster_key, cluster);

        std::cout << "Checking  si Track " << _track_map_silicon->find(si_track) << " cluster " << cluster_key
                  << " in layer " << layer << " position " << global(0) << "  " << global(1) << "  " << global(2)
                  << " eta " << si_track->get_eta() << std::endl;
      }
      else
      {
        std::cout << "Checking  si Track " << _track_map_silicon->find(si_track) << " cluster " << cluster_key
                  << " in layer " << layer << " with eta " << si_track->get_eta() << std::endl;
      }
    }

    if (trkrid == TrkrDefs::inttId)
    {
      TrkrCluster *cluster = _cluster_map->findCluster(cluster_key);
      if (!cluster)
      {
        continue;
      }

      unsigned int layer = TrkrDefs::getLayer(cluster_key);

      // get the bunch crossings for all hits in this cluster
      auto crossings = _cluster_crossing_map->getCrossings(cluster_key);
      for (auto iter1 = crossings.first; iter1 != crossings.second; ++iter1)
      {
        if (Verbosity() > 1)
        {
          std::cout << "                si Track " << _track_map_silicon->find(si_track) << " cluster " << iter1->first << " layer " << layer << " crossing " << iter1->second << std::endl;
        }
        intt_crossings.push_back(iter1->second);
      }
    }
  }

  return intt_crossings;
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
        std::cout << " drop si_track " << si_id << " with eta " << si_track->get_eta() << " and z " << si_track->get_z() << " because crossing is undefined " << std::endl;
      }
      continue;
    }

    float z_si = si_track->get_z();
    float z_tpc = tpc_track->get_z();
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

  //  std::cout << " tpc_pt = " << tpc_pt << " mag " << mag << " a " << match_function_a << " b " << match_function_b << std::endl;

  return mag;
}
