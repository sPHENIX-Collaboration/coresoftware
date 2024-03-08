#include "FillTruthRecoMatchTree.h"

#include <TFile.h>
#include <TH2D.h>
#include <TTree.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4tracking/EmbRecoMatch.h>
#include <g4tracking/EmbRecoMatchContainer.h>
#include <g4tracking/TrkrTruthTrack.h>
#include <g4tracking/TrkrTruthTrackContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxPHG4ParticleMap_v1.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <iostream>

using std::cout;
using std::endl;
//____________________________________________________________________________..
FillTruthRecoMatchTree::FillTruthRecoMatchTree(
    bool _fill_clusters, bool _fill_SvUnMatched, float _cluster_nzwidths, float _cluster_nphiwidths, const std::string& _outfile_name)
  : m_cluster_comp{_cluster_nphiwidths, _cluster_nzwidths}
  , m_fill_clusters{_fill_clusters}
  , m_fill_SvU{_fill_SvUnMatched}
  , m_outfile_name{_outfile_name}
{
  m_cluscntr.set_comparer(&m_cluster_comp);
  PHTFileServer::get().open(m_outfile_name, "RECREATE");

  h2_G4_nPixelsPhi = new TH2D("G4_nPixelsPhi", "PHG4 Emb Tracks; cluster pixel width Phi; layer",
                              100, -0.5, 99.5, 56, -0.5, 55.5);
  h2_G4_nPixelsZ = new TH2D("G4_nPixelsZ", "PHG4 Emb Tracks; cluster pixel width Z; layer",
                            100, -0.5, 99.5, 56, -0.5, 55.5);
  h2_Sv_nPixelsPhi = new TH2D("Sv_nPixelsPhi", "Svtx Reco Tracks; cluster pixel width Phi; layer",
                              100, -0.5, 99.5, 56, -0.5, 55.5);
  h2_Sv_nPixelsZ = new TH2D("Sv_nPixelsZ", "Svtx Reco Tracks; cluster pixel width Z; layer",
                            100, -0.5, 99.5, 56, -0.5, 55.5);

  m_ttree = new TTree("T", "Tracks (and sometimes clusters)");

  m_ttree->Branch("event", &nevent);
  m_ttree->Branch("nphg4_part", &nphg4_part);
  m_ttree->Branch("centrality", &centrality);
  m_ttree->Branch("ntrackmatches", &ntrackmatches);
  m_ttree->Branch("nphg4", &nphg4);
  m_ttree->Branch("nsvtx", &nsvtx);

  m_ttree->Branch("trackid", &b_trackid);
  m_ttree->Branch("is_G4track", &b_is_g4track);
  m_ttree->Branch("is_Svtrack", &b_is_Svtrack);
  m_ttree->Branch("is_matched", &b_is_matched);

  m_ttree->Branch("trkpt", &b_trkpt);
  m_ttree->Branch("trketa", &b_trketa);
  m_ttree->Branch("trkphi", &b_trkphi);

  m_ttree->Branch("nclus", &b_nclus);
  m_ttree->Branch("nclustpc", &b_nclustpc);
  m_ttree->Branch("nclusmvtx", &b_nclusmvtx);
  m_ttree->Branch("nclusintt", &b_nclusintt);
  m_ttree->Branch("matchrat", &b_matchrat);
  m_ttree->Branch("matchrat_intt", &b_matchrat_intt);
  m_ttree->Branch("matchrat_mvtx", &b_matchrat_mvtx);
  m_ttree->Branch("matchrat_tpc", &b_matchrat_tpc);

  if (m_fill_clusters)
  {
    m_ttree->Branch("clus_match", &b_clusmatch);
    m_ttree->Branch("clus_x", &b_clus_x);
    m_ttree->Branch("clus_y", &b_clus_y);
    m_ttree->Branch("clus_z", &b_clus_z);
    m_ttree->Branch("clus_r", &b_clus_r);
    m_ttree->Branch("clus_layer", &b_clus_layer);
    m_ttree->Branch("nphibins", &b_clus_nphibins);
    m_ttree->Branch("nzbins", &b_clus_ntbins);
  }
}

//____________________________________________________________________________..
FillTruthRecoMatchTree::~FillTruthRecoMatchTree() = default;

//____________________________________________________________________________..
int FillTruthRecoMatchTree::Init(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << " Beginning FillTruthRecoMatchTree " << std::endl;
    topNode->print();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int FillTruthRecoMatchTree::InitRun(PHCompositeNode* topNode)
{
  auto init_status = m_cluster_comp.init(topNode);
  if (init_status == Fun4AllReturnCodes::ABORTRUN)
  {
    return init_status;
  }

  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int FillTruthRecoMatchTree::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cout << PHWHERE << " DST node is missing, quitting" << std::endl;
    std::cerr << PHWHERE << " DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in FillTruthRecoMatchTree::createNodes");
  }

  m_EmbRecoMatchContainer = findNode::getClass<EmbRecoMatchContainer>(topNode, "TRKR_EMBRECOMATCHCONTAINER");
  if (!m_EmbRecoMatchContainer)
  {
    std::cout << PHWHERE << " Cannot find node TRKR_EMBRECOMATCHCONTAINER on node tree; quitting " << std::endl;
    std::cerr << PHWHERE << " Cannot find node TRKR_EMBRECOMATCHCONTAINER on node tree; quitting " << std::endl;
    throw std::runtime_error(" Cannot find node TRKR_EMBRECOMATCHCONTAINER on node tree; quitting");
  }

  PHCompositeNode* svtxNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_PHG4TruthInfoContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_PHG4TruthInfoContainer)
  {
    std::cout << "Could not locate G4TruthInfo node when running "
              << "\"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_SvtxTrackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!m_SvtxTrackMap)
  {
    std::cout << "Could not locate SvtxTrackMap node when running "
              << "\"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  /* m_TruthClusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, */
  /*     "TRKR_TRUTHCLUSTERCONTAINER"); */
  /* if (!m_TruthClusterContainer) */
  /* { */
  /*   std::cout << "Could not locate TRKR_TRUTHCLUSTERCONTAINER node when running " */
  /*     << "\"TruthRecoTrackMatching\" module." << std::endl; */
  /*   return Fun4AllReturnCodes::ABORTEVENT; */
  /* } */

  /* m_RecoClusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER"); */
  /* if (!m_RecoClusterContainer) */
  /* { */
  /*   std::cout << "Could not locate TRKR_CLUSTER node when running " */
  /*     << "\"TruthRecoTrackMatching\" module." << std::endl; */
  /*   return Fun4AllReturnCodes::ABORTEVENT; */
  /* } */

  m_TrkrTruthTrackContainer = findNode::getClass<TrkrTruthTrackContainer>(topNode,
                                                                          "TRKR_TRUTHTRACKCONTAINER");
  if (!m_TrkrTruthTrackContainer)
  {
    std::cout << "Could not locate TRKR_TRUTHTRACKCONTAINER node when running "
              << "\"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int FillTruthRecoMatchTree::process_event(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 5)
  {
    cout << " FillTruthRecoMatchTree::process_event() " << endl;
  }

  // fill in the event data
  ++nevent;
  nphg4 = m_TrkrTruthTrackContainer->getMap().size();
  nsvtx = m_SvtxTrackMap->size();
  ntrackmatches = m_EmbRecoMatchContainer->getMatches().size();
  // get centrality later...

  // fill in pixel widths on truth tracks
  for (auto hitsetkey : m_cluscntr.get_PHG4_clusters()->getHitSetKeys())
  {
    float layer = (float) TrkrDefs::getLayer(hitsetkey);
    auto range = m_cluscntr.get_PHG4_clusters()->getClusters(hitsetkey);
    for (auto iter = range.first; iter != range.second; ++iter)
    {
      auto& cluster = iter->second;
      h2_G4_nPixelsPhi->Fill((float) cluster->getPhiSize(), layer);
      h2_G4_nPixelsZ->Fill((float) cluster->getZSize(), layer);
    }
  }
  // fill in pixel widths on reco tracks
  for (auto hitsetkey : m_cluscntr.get_SVTX_clusters()->getHitSetKeys())
  {
    float layer = (float) TrkrDefs::getLayer(hitsetkey);
    auto range = m_cluscntr.get_SVTX_clusters()->getClusters(hitsetkey);
    for (auto iter = range.first; iter != range.second; ++iter)
    {
      auto& cluster = iter->second;
      h2_Sv_nPixelsPhi->Fill((float) cluster->getPhiSize(), layer);
      h2_Sv_nPixelsZ->Fill((float) cluster->getZSize(), layer);
    }
  }

  nphg4_part = 0;
  const auto range = m_PHG4TruthInfoContainer->GetPrimaryParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second; ++iter)
  {
    nphg4_part++;
  }

  // unmatched tracks are only entered once
  // matches can repeat a given svtx or phg4 track, depending on the
  // parameters in teh matching in filltruthrecomatchtree
  //
  // (1) fill unmatched phg4
  // (2) fill unmatched svtx
  // (3) fill matched phg4 and svtx
  clear_clusvecs(" nothing ");

  if (Verbosity() > 2)
  {
    std::cout << " getting" << (int) m_EmbRecoMatchContainer->getMatches().size() << std::endl;
  }
  for (auto& match : m_EmbRecoMatchContainer->getMatches())
  {
    unsigned int g4_trkid = match->idTruthTrack();
    int sv_trkid = match->idRecoTrack();

    auto g4trk = m_TrkrTruthTrackContainer->getTruthTrack(g4_trkid);
    auto svtrk = m_SvtxTrackMap->get(sv_trkid);

    m_cluscntr.addClusKeys(g4trk);
    m_cluscntr.addClusKeys(svtrk);
    m_cluscntr.find_matches();

    // <- <- <- <- G4 Matched Tracks
    b_is_matched = true;
    b_is_g4track = true;
    b_is_Svtrack = false;

    b_trackid = g4_trkid;
    b_trkpt = g4trk->getPt();
    b_trketa = g4trk->getPseudoRapidity();
    b_trkphi = g4trk->getPhi();

    auto cnt = m_cluscntr.phg4_cntclus();
    auto cnt_match = m_cluscntr.phg4_cnt_matchedclus();

    b_nclus = cnt[4];
    b_nclusmvtx = cnt[0];
    b_nclusintt = cnt[1];
    b_nclustpc = cnt[2];

    b_matchrat = (float) cnt_match[4] / cnt[4];
    b_matchrat_mvtx = (float) cnt_match[0] / cnt[0];
    b_matchrat_intt = (float) cnt_match[1] / cnt[1];
    b_matchrat_tpc = (float) cnt_match[2] / cnt[2];

    if (m_fill_clusters)
    {
      auto clusters = m_cluscntr.phg4_clusloc_unmatched();
      for (auto& loc : clusters)
      {
        b_clusmatch.push_back(false);
        b_clus_layer.push_back(std::get<0>(loc));
        auto x = std::get<1>(loc)[0];
        auto y = std::get<1>(loc)[1];
        b_clus_x.push_back(x);
        b_clus_y.push_back(y);
        b_clus_z.push_back(std::get<1>(loc)[2]);
        b_clus_r.push_back(pow(x * x + y * y, 0.5));
        b_clus_nphibins.push_back(std::get<2>(loc));
        b_clus_ntbins.push_back(std::get<3>(loc));
      }

      clusters = m_cluscntr.clusloc_matched();
      for (auto& loc : clusters)
      {
        b_clusmatch.push_back(true);
        b_clus_layer.push_back(std::get<0>(loc));
        auto x = std::get<1>(loc)[0];
        auto y = std::get<1>(loc)[1];
        b_clus_x.push_back(x);
        b_clus_y.push_back(y);
        b_clus_z.push_back(std::get<1>(loc)[2]);
        b_clus_r.push_back(pow(x * x + y * y, 0.5));
        b_clus_nphibins.push_back(std::get<2>(loc));
        b_clus_ntbins.push_back(std::get<3>(loc));
      }
    }
    m_ttree->Fill();
    clear_clusvecs();
    /* clear_clusvecs("apple0 g4_matched"); */

    // <- <- <- <- Svtx Matched Tracks
    b_is_g4track = false;
    b_is_Svtrack = true;
    b_trackid = sv_trkid;
    b_trkpt = svtrk->get_pt();
    b_trketa = svtrk->get_eta();
    b_trkphi = svtrk->get_phi();

    cnt = m_cluscntr.svtx_cntclus();
    b_nclus = cnt[4];
    b_nclusmvtx = cnt[0];
    b_nclusintt = cnt[1];
    b_nclustpc = cnt[2];

    b_matchrat = (float) cnt_match[4] / cnt[4];
    b_matchrat_mvtx = (float) cnt_match[0] / cnt[0];
    b_matchrat_intt = (float) cnt_match[1] / cnt[1];
    b_matchrat_tpc = (float) cnt_match[2] / cnt[2];

    /* int _ = 0; */
    if (m_fill_clusters)
    {
      auto clusters = m_cluscntr.svtx_clusloc_unmatched();
      for (auto& loc : clusters)
      {
        b_clusmatch.push_back(false);
        b_clus_layer.push_back(std::get<0>(loc));
        auto x = std::get<1>(loc)[0];
        auto y = std::get<1>(loc)[1];
        /* if (_==0) cout << " apple x: " << x << " y: " << y << endl; */
        /* _ += 1; */
        b_clus_x.push_back(x);
        b_clus_y.push_back(y);
        b_clus_z.push_back(std::get<1>(loc)[2]);
        b_clus_r.push_back(pow(x * x + y * y, 0.5));
        b_clus_nphibins.push_back(std::get<2>(loc));
        b_clus_ntbins.push_back(std::get<3>(loc));
      }

      clusters = m_cluscntr.clusloc_matched();
      for (auto& loc : clusters)
      {
        b_clusmatch.push_back(true);
        b_clus_layer.push_back(std::get<0>(loc));
        auto x = std::get<1>(loc)[0];
        auto y = std::get<1>(loc)[1];
        b_clus_x.push_back(x);
        b_clus_y.push_back(y);
        b_clus_z.push_back(std::get<1>(loc)[2]);
        b_clus_r.push_back(pow(x * x + y * y, 0.5));
        b_clus_nphibins.push_back(std::get<2>(loc));
        b_clus_ntbins.push_back(std::get<3>(loc));
      }
    }
    m_ttree->Fill();
    clear_clusvecs();
    /* clear_clusvecs("apple1 s4_matched"); */
  }

  // <- <- <- <- G4 un-matched Tracks
  b_is_matched = false;
  b_is_g4track = true;
  b_is_Svtrack = false;
  for (auto& g4_trkid : m_EmbRecoMatchContainer->ids_TruthUnmatched())
  {
    auto g4trk = m_TrkrTruthTrackContainer->getTruthTrack(g4_trkid);
    m_cluscntr.addClusKeys(g4trk);

    b_trackid = g4_trkid;
    b_trkpt = g4trk->getPt();
    b_trketa = g4trk->getPseudoRapidity();
    b_trkphi = g4trk->getPhi();

    auto cnt = m_cluscntr.phg4_cntclus();
    b_nclus = cnt[4];
    b_nclusmvtx = cnt[0];
    b_nclusintt = cnt[1];
    b_nclustpc = cnt[2];

    if (m_fill_clusters)
    {
      auto clusters = m_cluscntr.phg4_clusloc_all();
      for (auto& loc : clusters)
      {
        b_clusmatch.push_back(false);
        b_clus_layer.push_back(std::get<0>(loc));
        auto x = std::get<1>(loc)[0];
        auto y = std::get<1>(loc)[1];
        b_clus_x.push_back(x);
        b_clus_y.push_back(y);
        b_clus_z.push_back(std::get<1>(loc)[2]);
        b_clus_r.push_back(pow(x * x + y * y, 0.5));
        b_clus_nphibins.push_back(std::get<2>(loc));
        b_clus_ntbins.push_back(std::get<3>(loc));
      }
      // this is an unmatched track, so there are no matched clusters
    }
    m_ttree->Fill();
    clear_clusvecs();
    /* clear_clusvecs("apple2 g4_unmatched"); */
  }

  // <- <- <- <- Svtx un-matched Tracks
  b_is_matched = false;
  b_is_matched = false;
  b_is_g4track = false;
  b_is_Svtrack = true;

  // just put in all svtx tracks, period...
  if (m_fill_SvU)
  {
    for (auto sv_trkid : G4Eval::unmatchedSvtxTrkIds(m_EmbRecoMatchContainer, m_SvtxTrackMap))
    {
      auto svtrk = m_SvtxTrackMap->get(sv_trkid);
      m_cluscntr.addClusKeys(svtrk);
      b_trackid = sv_trkid;
      b_trkpt = svtrk->get_pt();
      b_trketa = svtrk->get_eta();
      b_trkphi = svtrk->get_phi();

      auto cnt = m_cluscntr.svtx_cntclus();
      b_nclus = cnt[4];
      b_nclusmvtx = cnt[0];
      b_nclusintt = cnt[1];
      b_nclustpc = cnt[2];

      if (m_fill_clusters)
      {
        auto clusters = m_cluscntr.svtx_clusloc_all();
        for (auto& loc : clusters)
        {
          b_clusmatch.push_back(false);
          b_clus_layer.push_back(std::get<0>(loc));
          auto x = std::get<1>(loc)[0];
          auto y = std::get<1>(loc)[1];
          b_clus_x.push_back(x);
          b_clus_y.push_back(y);
          b_clus_z.push_back(std::get<1>(loc)[2]);
          b_clus_r.push_back(pow(x * x + y * y, 0.5));
          b_clus_nphibins.push_back(std::get<2>(loc));
          b_clus_ntbins.push_back(std::get<3>(loc));
        }
      }
      m_ttree->Fill();
      clear_clusvecs();
    }
  }

  if (Verbosity() > 100)
  {
    print_mvtx_diagnostics();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void FillTruthRecoMatchTree::print_mvtx_diagnostics()
{
  std::cout << "To do: "
            << " (1)  number of truth tracks and total number of mvtx and ratio " << std::endl
            << " (2)  ditto for reco tracks " << std::endl;

  double n_PHG4_tracks = m_TrkrTruthTrackContainer->getMap().size();
  // count how many mvtx clusters in the phg4 tracks
  double n_in_PHG4_tracks{0.};
  for (auto& pair : m_TrkrTruthTrackContainer->getMap())
  {
    m_cluscntr.addClusKeys(pair.second);
    n_in_PHG4_tracks += m_cluscntr.phg4_cntclus()[0];
  }
  // count how mant mvtx clusters in truth container (should be the same)
  double n_in_PHG4_clusters{0.};
  for (auto hitsetkey : m_cluscntr.get_PHG4_clusters()->getHitSetKeys())
  {
    if (TrkrDefs::getLayer(hitsetkey) > 2)
    {
      continue;
    }
    auto range = m_cluscntr.get_PHG4_clusters()->getClusters(hitsetkey);
    for (auto r = range.first; r != range.second; ++r)
    {
      n_in_PHG4_clusters += 1.;
    }
  }

  // count how many svtx tracks
  double n_SVTX_tracks = m_SvtxTrackMap->size();
  // count how many mvtx clusters in svtx tracks
  double n_in_SVTX_tracks{0.};
  for (auto& entry : *m_SvtxTrackMap)
  {
    m_cluscntr.addClusKeys(entry.second);
    n_in_SVTX_tracks += m_cluscntr.svtx_cntclus()[0];
  }
  // count how many mvtx are total in the container
  double n_in_SVTX_clusters{0.};
  for (auto hitsetkey : m_cluscntr.get_SVTX_clusters()->getHitSetKeys())
  {
    if (TrkrDefs::getLayer(hitsetkey) > 2)
    {
      continue;
    }
    auto range = m_cluscntr.get_SVTX_clusters()->getClusters(hitsetkey);
    for (auto r = range.first; r != range.second; ++r)
    {
      n_in_SVTX_clusters += 1.;
    }
  }

  std::cout << Form(
                   "MVTX"
                   "\nPHG4:  Tracks(%.0f)   Clusters In tracks(%.0f)   Total (%.0f)"
                   "\n       ave. per track: %6.3f   ratio in all tracks: %6.2f",
                   n_PHG4_tracks, n_in_PHG4_tracks, n_in_PHG4_clusters, (n_in_PHG4_tracks / n_PHG4_tracks), (n_in_PHG4_tracks / n_in_PHG4_clusters))
            << std::endl;
  std::cout << Form(
                   "\nSVTX:  Tracks(%.0f)   Clusters In tracks(%.0f)   Total (%.0f)"
                   "\n       ave. per track: %6.3f   ratio in all tracks: %6.2f",
                   n_SVTX_tracks, n_in_SVTX_tracks, n_in_SVTX_clusters, (n_in_SVTX_tracks / n_SVTX_tracks), (n_in_SVTX_tracks / n_in_SVTX_clusters))
            << std::endl;
}

int FillTruthRecoMatchTree::End(PHCompositeNode* /*unused*/)
{
  if (Verbosity() > 2)
  {
    std::cout << PHWHERE << ": ending FillTruthRecoMatchTree" << std::endl;
  }
  PHTFileServer::get().cd(m_outfile_name);

  h2_G4_nPixelsPhi->Write();
  h2_G4_nPixelsZ->Write();
  h2_Sv_nPixelsPhi->Write();
  h2_Sv_nPixelsZ->Write();

  m_ttree->Write();
  return Fun4AllReturnCodes::EVENT_OK;
}

void FillTruthRecoMatchTree::clear_clusvecs(const std::string& tag)
{
  /* cout << " banana |" << tag << "|"<<endl; */
  if (tag != "")
  {
    cout << endl
         << " pear printing " << tag << " x(" << b_clus_x.size() << ") ";
    for (auto x : b_clus_x)
    {
      cout << x << " ";
    }
    cout << endl;
  }
  // Tracks and clustes
  b_clusmatch.clear();
  b_clus_x.clear();
  b_clus_y.clear();
  b_clus_z.clear();
  b_clus_r.clear();
  b_clus_layer.clear();
  b_clus_nphibins.clear();
  b_clus_ntbins.clear();
}
