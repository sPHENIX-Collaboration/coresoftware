#include "MvtxMatchingEfficiencyWithShapes.h"

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackAnalysisUtils.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/SvtxTrack.h>              // for SvtxTrack
#include <trackbase_historic/SvtxTrackState.h>         // for SvtxTrackState

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TH1.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>                                    // for basic_ostream
#include <set>
#include <utility>                                     // for get, pair

namespace
{
  //! range adaptor to be able to use range-based for loop
  template <class T>
  class range_adaptor
  {
   public:
    explicit range_adaptor(const T& range)
      : m_range(range)
    {
    }
    const typename T::first_type& begin() { return m_range.first; }
    const typename T::second_type& end() { return m_range.second; }

   private:
    T m_range;
  };
}  // namespace

//____________________________________________________________________________..
MvtxMatchingEfficiencyWithShapes::MvtxMatchingEfficiencyWithShapes(const std::string& name, const std::string& outputfilename)
  : SubsysReco(name)
  , m_outputFileName(outputfilename)
{
}

//____________________________________________________________________________..
int MvtxMatchingEfficiencyWithShapes::InitRun(PHCompositeNode* /*topNode*/)
{
  // print configuration
  std::cout << "MvtxMatchingEfficiencyWithShapes::InitRun " << std::endl;

  PHTFileServer::open(m_outputFileName, "RECREATE");

  h_status = new TH1I("h_status_wTPC", "h_status_wTPC", 8, 0, 8);
  h_status->GetXaxis()->SetBinLabel(1, "Total");
  h_status->GetXaxis()->SetBinLabel(2, "INTT >= 2");
  h_status->GetXaxis()->SetBinLabel(3, "INTT same crossing");
  h_status->GetXaxis()->SetBinLabel(4, "INTT unique clusters");
  h_status->GetXaxis()->SetBinLabel(5, "MVTX unique clusters");
  h_status->GetXaxis()->SetBinLabel(6, "z < 10 cm");
  h_status->GetXaxis()->SetBinLabel(7, "eta < 1.1");
  h_status->GetXaxis()->SetBinLabel(8, "MVTX >= 2 layer");

  h_INTT_time_delta = new TH1I("h_INTT_time_delta_wTPC", "h_INTT_time_delta_wTPC", 900, 0.5, 900.5);

  // Create a new TTree
  tree = new TTree("TPC_tracks", "TPC_tracks");

  // Create branches
  tree->Branch("pt", &pt, "pt/F");
  tree->Branch("eta", &eta, "eta/F");
  tree->Branch("phi", &phi, "phi/F");
  tree->Branch("frac_p_z", &frac_p_z, "frac_p_z/F");
  tree->Branch("dEdx", &dEdx, "dEdx/F");
  tree->Branch("nTPC", &nTPC, "nTPC/I");
  tree->Branch("layers", &layers, "layers/I");
  tree->Branch("states", &states, "states/I");
// NOLINTNEXTLINE (readability-container-data-pointer)
  tree->Branch("shape_L0_C0_nhits", nhits_arr[0].data(), "shape_L0_C0_nhits/I");
  tree->Branch("shape_L0_C0_x", x_arr[0].data(), "shape_L0_C0_x/I");
  tree->Branch("shape_L0_C0_y", y_arr[0].data(), "shape_L0_C0_y/I");
  tree->Branch("shape_L0_C0_key", key_arr[0].data(), "shape_L0_C0_key/l");

  tree->Branch("shape_L0_C1_nhits", &nhits_arr[0][1], "shape_L0_C1_nhits/I");
  tree->Branch("shape_L0_C1_x", &x_arr[0][1], "shape_L0_C1_x/I");
  tree->Branch("shape_L0_C1_y", &y_arr[0][1], "shape_L0_C1_y/I");
  tree->Branch("shape_L0_C1_key", &key_arr[0][1], "shape_L0_C1_key/l");

  tree->Branch("shape_L1_C0_nhits", nhits_arr[1].data(), "shape_L1_C0_nhits/I");
  tree->Branch("shape_L1_C0_x", x_arr[1].data(), "shape_L1_C0_x/I");
  tree->Branch("shape_L1_C0_y", y_arr[1].data(), "shape_L1_C0_y/I");
  tree->Branch("shape_L1_C0_key", key_arr[1].data(), "shape_L1_C0_key/l");

  tree->Branch("shape_L1_C1_nhits", &nhits_arr[1][1], "shape_L1_C1_nhits/I");
  tree->Branch("shape_L1_C1_x", &x_arr[1][1], "shape_L1_C1_x/I");
  tree->Branch("shape_L1_C1_y", &y_arr[1][1], "shape_L1_C1_y/I");
  tree->Branch("shape_L1_C1_key", &key_arr[1][1], "shape_L1_C1_key/l");

  tree->Branch("shape_L2_C0_nhits", nhits_arr[2].data(), "shape_L2_C0_nhits/I");
  tree->Branch("shape_L2_C0_x", x_arr[2].data(), "shape_L2_C0_x/I");
  tree->Branch("shape_L2_C0_y", y_arr[2].data(), "shape_L2_C0_y/I");
  tree->Branch("shape_L2_C0_key", key_arr[2].data(), "shape_L2_C0_key/l");

  tree->Branch("shape_L2_C2_nhits", &nhits_arr[2][1], "shape_L2_C2_nhits/I");
  tree->Branch("shape_L2_C2_x", &x_arr[2][1], "shape_L2_C2_x/I");
  tree->Branch("shape_L2_C2_y", &y_arr[2][1], "shape_L2_C2_y/I");
  tree->Branch("shape_L2_C2_key", &key_arr[2][1], "shape_L2_C2_key/l");

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MvtxMatchingEfficiencyWithShapes::process_event(PHCompositeNode* topNode)
{
  ievent++;

  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!trackmap)
  {
    std::cout << PHWHERE
              << "SvtxTrackMap node is missing, can't collect particles"
              << std::endl;
    return -1;
  }

  TrackSeedContainer* silicon_track_map = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if (!silicon_track_map)
  {
    std::cout << PHWHERE
              << "SiliconTrackSeedContainer node is missing, can't collect particles"
              << std::endl;
    return -1;
  }

  m_hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hitsetcontainer)
  {
    std::cout << "m_hitsetcontainer not found" << std::endl;
  }
  // assert(m_hitsetcontainer);

  // cluster hit association map
  m_cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!m_cluster_hit_map)
  {
    std::cout << "m_cluster_hit_map not found" << std::endl;
  }

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
  if (!_cluster_map)
  {
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  }

  if (!_cluster_map)
  {
    std::cout << PHWHERE
              << "TrkrClusterContainer node is missing"
              << std::endl;
    return -1;
  }

  _geom_container = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  if (!_geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node TPCGEOMCONTAINER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE << "No acts reco geometry, bailing.";
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // std::cout<<"SvtxTrackMap"<<std::endl;
  // fetch MVTX clusters
  std::vector<TrkrDefs::cluskey> cluster_vector;
  std::vector<TrkrDefs::cluskey> cluster_vector_INTT;
  for (auto& it : *trackmap)
  {
    for (const auto& ckey : TrackAnalysisUtils::get_cluster_keys(it.second))
    {
      switch (TrkrDefs::getTrkrId(ckey))
      {
      case TrkrDefs::mvtxId:
        cluster_vector.push_back(ckey);
        break;
      case TrkrDefs::inttId:
        cluster_vector_INTT.push_back(ckey);
        break;
      default:
	break;
      }
    }
  }

  // find duplicated MVTX clusters in an event
  std::set<TrkrDefs::cluskey> duplicate_cluster_vector = findDuplicates(cluster_vector);
  // find duplicated INTT clusters in an event
  std::set<TrkrDefs::cluskey> duplicate_cluster_vector_INTT = findDuplicates(cluster_vector_INTT);

  for (auto& it : *trackmap)
  {
    // number of cluster per subsystem
    int n_intt_hits = 0;
    int n_tpc_hits = 0;
    // fired mvtx layers 0 - innermost , 2 - outermost
    bool mvtx_l[3] = {false, false, false};
    bool mvtx_l_state[3] = {false, false, false};
    std::vector<float> intt_time;

    SvtxTrack* track = it.second;
    // count clusters per subsystem and fired mvtx layers
    for (const auto& ckey : TrackAnalysisUtils::get_cluster_keys(track))
    {
      switch (TrkrDefs::getTrkrId(ckey))
      {
      case TrkrDefs::mvtxId:
        mvtx_l[TrkrDefs::getLayer(ckey)] = true;
        break;
      case TrkrDefs::inttId:
        intt_time.push_back(InttDefs::getTimeBucketId(ckey));
        n_intt_hits++;
        break;
      case TrkrDefs::tpcId:
        n_tpc_hits++;
        break;
      default:
	break;
      }
    }

    // check states in MVTX layers
    for (auto state_it = track->begin_states(); state_it != track->end_states(); ++state_it)
    {
      auto clus_key = state_it->second->get_cluskey();
      switch (TrkrDefs::getTrkrId(clus_key))
      {
      case TrkrDefs::mvtxId:
        mvtx_l_state[TrkrDefs::getLayer(clus_key)] = true;
        break;
      default:
	break;
      }
    }

    // all events
    h_status->Fill(0.5);

    // require number of INTT clusters >= 2
    if (n_intt_hits < 2)
    {
      continue;
    }
    h_status->Fill(1.5);

    // require the INTT clusters are from same crossing
    bool good_intt_time = std::all_of(intt_time.begin() + 1, intt_time.end(), [&](int i)
                                      { return i == intt_time[0]; });

    if (good_intt_time == false)
    {
      if (n_intt_hits == 2)
      {
        h_INTT_time_delta->Fill(std::abs(intt_time[0] - intt_time[1]));
      }
      continue;
    }
    intt_time.clear();
    h_status->Fill(2.5);

    // require the tracks does not share INTT clusters
    bool good_track_INTT = true;
    for (const auto& ckey : TrackAnalysisUtils::get_cluster_keys(track))
    {
      if (duplicate_cluster_vector_INTT.contains(ckey))
      {
        good_track_INTT = false;
      }
    }
    if (good_track_INTT == false)
    {
      continue;
    }
    h_status->Fill(3.5);

    // require the tracks does not share MVTX clusters
    bool good_track = true;
    for (const auto& ckey : TrackAnalysisUtils::get_cluster_keys(track))
    {
      if (duplicate_cluster_vector.contains(ckey))
      {
        good_track = false;
      }
    }
    if (good_track == false)
    {
      continue;
    }
    h_status->Fill(4.5);

    // z position cut
    if (std::abs(track->get_z()) > 10)
    {
      continue;
    }
    h_status->Fill(5.5);

    // eta cut
    if (std::abs(track->get_eta()) > 1.1)
    {
      continue;
    }
    h_status->Fill(6.5);

    // require 2 or 3 MVTX layers fired (layers not clusters)
    if (std::count_if(std::begin(mvtx_l), std::end(mvtx_l), [](bool i)
                      { return i == true; }) < 2)
    {
      continue;
    }
    h_status->Fill(7.5);

    // fill Ttree variables
    pt = track->get_pt();
    eta = track->get_eta();
    phi = track->get_phi();
    frac_p_z = track->get_p() / track->get_charge();
    dEdx = calc_dedx(track->get_tpc_seed());
    layers = -1;
    states = -1;
    nTPC = n_tpc_hits;

    // initiate cluster shape variables
    int layer_cluster_counter[3] = {-1, -1, -1};
    for (auto& inner_array : nhits_arr)
    {
      std::fill(inner_array.begin(), inner_array.end(), 0);
    }
    for (auto& inner_array : x_arr)
    {
      std::fill(inner_array.begin(), inner_array.end(), 0);
    }
    for (auto& inner_array : y_arr)
    {
      std::fill(inner_array.begin(), inner_array.end(), 0);
    }
    for (auto& inner_array : key_arr)
    {
      std::fill(inner_array.begin(), inner_array.end(), 0);
    }

    // calculate cluster shape
    for (const auto& ckey : TrackAnalysisUtils::get_cluster_keys(track))
    {
      switch (TrkrDefs::getTrkrId(ckey))
      {
      case TrkrDefs::mvtxId:
      {
        int layer = TrkrDefs::getLayer(ckey);
        // std::cout << "layer " << layer << std::endl;
        layer_cluster_counter[layer]++;
        int cluster_index = layer_cluster_counter[layer];

        // find associated hits
        const auto hit_range = m_cluster_hit_map->getHits(ckey);

        int cL = 2000;
        int cR = 0;
        int cT = 2000;
        int cB = 0;
        uint64_t cmap = 0;

        for (const auto& [ckey2, hitkey] : range_adaptor(hit_range))
        {
          int col = MvtxDefs::getCol(hitkey);
          int row = MvtxDefs::getRow(hitkey);
          cL = std::min(cL, col);
          cR = std::max(cR, col);
          cT = std::min(cT, row);
          cB = std::max(cB, row);
        }

        int x = cR - cL + 1;
        int y = cB - cT + 1;
        int nhits = std::distance(hit_range.first, hit_range.second);

        // key is saved only for shapes < 64 hits
        if (x * y < 64)
        {
          for (const auto& [ckey2, hitkey] : range_adaptor(hit_range))
          {
	    // NOLINTNEXTLINE(hicpp-signed-bitwise)
            cmap |= (1ULL << ((MvtxDefs::getCol(hitkey) - cL) + ((cR - cL + 1) * (MvtxDefs::getRow(hitkey) - cT))));
            // std::cout<<"toggle bit "<<((MvtxDefs::getCol(hitkey)-cL)+((cR-cL+1)*(MvtxDefs::getRow(hitkey)-cT)))<<" "<< (MvtxDefs::getCol(hitkey)-cL)<<" "<<(cR-cL+1) <<" "<< (-MvtxDefs::getRow(hitkey)+cT)<<std::endl;
          }
          // std::cout<<"key "<<cmap<<std::endl;
        }

        // shape debug prints

        // std::cout<<"L: "<<cL<<" R: "<<cR<<" B: "<<cB<<" T: "<<cT<<std::endl;
        // std::cout<<"key int: "<<cmap<<std::endl;
        // std::cout<<"key uint: "<<(unsigned long long)cmap<<std::endl;
        // std::cout<<"key int: ";

        // std::bitset<64> bits(cmap);
        // std::cout << bits << std::endl;

        // Using bitwise operations
        // for (int i = 63; i >= 0; --i) {
        //    std::cout << ((cmap >> i) & 1);
        //}
        // std::cout << std::endl;

        // for( const auto& [ckey2, hitkey]:range_adaptor( hit_range ) ){
        //   std::cout<<"row: "<<MvtxDefs::getRow(hitkey)<<" col "<< MvtxDefs::getCol(hitkey)<<" bitmap "<<((MvtxDefs::getCol(hitkey)-cL)+((cR-cL+1)*(MvtxDefs::getRow(hitkey)-cT)))<<std::endl;
        // }

        // Fill the arrays
        nhits_arr[layer][cluster_index] = nhits;
        x_arr[layer][cluster_index] = x;
        y_arr[layer][cluster_index] = y;
        key_arr[layer][cluster_index] = cmap;

        break;
      }
      default:
	break;
      }
    }

    // encode layer map
    layers = (mvtx_l[0] ? 1U : 0)     // L0 → bit 0
             | (mvtx_l[1] ? 2U : 0)   // L1 → bit 1
             | (mvtx_l[2] ? 4U : 0);  // L2 → bit 2

    // encode state map
    states = (mvtx_l_state[0] ? 1U : 0)     // L0 → bit 0
             | (mvtx_l_state[1] ? 2U : 0)   // L1 → bit 1
             | (mvtx_l_state[2] ? 4U : 0);  // L2 → bit 2
    tree->Fill();
  }

  // clear event vectors for tracks
  cluster_vector.clear();
  duplicate_cluster_vector.clear();
  cluster_vector_INTT.clear();
  duplicate_cluster_vector_INTT.clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MvtxMatchingEfficiencyWithShapes::EndRun(const int /*runnumber*/)
{
  // TFile *root_out = new TFile("MVTX_ME.root","RECREATE");

  std::cout << "MvtxMatchingEfficiencyWithShapes::End - Output to " << m_outputFileName << std::endl;

  if (PHTFileServer::cd(m_outputFileName))
  {
    h_status->Write();
    h_INTT_time_delta->Write();
    tree->Write();
  }

  std::cout << "MvtxMatchingEfficiencyWithShapes::End(PHCompositeNode *topNode) This is the End..." << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

float MvtxMatchingEfficiencyWithShapes::calc_dedx(TrackSeed* tpcseed)
{
  std::vector<TrkrDefs::cluskey> clusterKeys;
  clusterKeys.insert(clusterKeys.end(), tpcseed->begin_cluster_keys(), tpcseed->end_cluster_keys());

  std::vector<float> dedxlist;
  for (unsigned long cluster_key : clusterKeys)
  {
    unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
    if (TrkrDefs::getTrkrId(cluster_key) != TrkrDefs::TrkrId::tpcId)
    {
      continue;
    }
    TrkrCluster* cluster = _cluster_map->findCluster(cluster_key);

    float adc = cluster->getAdc();
    PHG4TpcGeom* GeoLayer_local = _geom_container->GetLayerCellGeom(layer_local);
    float thick = GeoLayer_local->get_thickness();

    float r = GeoLayer_local->get_radius();
    float alpha = (r * r) / (2 * r * std::abs(1.0 / tpcseed->get_qOverR()));
    float beta = std::atan(tpcseed->get_slope());
    float alphacorr = std::cos(alpha);
    if (alphacorr < 0 || alphacorr > 4)
    {
      alphacorr = 4;
    }
    float betacorr = std::cos(beta);
    if (betacorr < 0 || betacorr > 4)
    {
      betacorr = 4;
    }
    adc /= thick;
    adc *= alphacorr;
    adc *= betacorr;
    dedxlist.push_back(adc);
    sort(dedxlist.begin(), dedxlist.end());
  }
  int trunc_min = 0;
  int trunc_max = (int) dedxlist.size() * 0.7;
  float sumdedx = 0;
  int ndedx = 0;
  for (int j = trunc_min; j <= trunc_max; j++)
  {
    sumdedx += dedxlist.at(j);
    ndedx++;
  }
  sumdedx /= ndedx;
  return sumdedx;
}

std::set<TrkrDefs::cluskey> MvtxMatchingEfficiencyWithShapes::findDuplicates(std::vector<TrkrDefs::cluskey> vec)
{
  std::set<TrkrDefs::cluskey> duplicates;
  std::sort(vec.begin(), vec.end());
  std::set<TrkrDefs::cluskey> distinct(vec.begin(), vec.end());
  std::set_difference(vec.begin(), vec.end(), distinct.begin(), distinct.end(),
                      std::inserter(duplicates, duplicates.end()));
  return duplicates;
}
