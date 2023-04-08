#include "FillTruthRecoMatchTree.h"

#include <TFile.h>
#include <TTree.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4tracking/EmbRecoMatch.h>
#include <g4tracking/EmbRecoMatchContainer.h>
#include <g4tracking/TrkrTruthTrack.h>
#include <g4tracking/TrkrTruthTrackContainer.h>
#include <iostream>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxPHG4ParticleMap_v1.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <TH2D.h>
#include <trackbase/TrkrCluster.h>

using std::cout;
using std::endl;
//____________________________________________________________________________..
FillTruthRecoMatchTree::FillTruthRecoMatchTree(
      bool _fill_clusters
    , bool _fill_SvUnMatched
    , float _cluster_nzwidths
    , float _cluster_nphiwidths
    , const std::string& _outfile_name
    )
  : 
    m_cluster_comp  { _cluster_nphiwidths, _cluster_nzwidths }
  , m_fill_clusters { _fill_clusters }
  , m_fill_SvU      { _fill_SvUnMatched }
  , m_outfile_name  { _outfile_name }
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

  m_ttree->Branch("event",                  &nevent);
  m_ttree->Branch("nphg4_part",             &nphg4_part);
  m_ttree->Branch("centrality",             &centrality);
  m_ttree->Branch("ntrackmatches",          &ntrackmatches);
  m_ttree->Branch("nphg4",                  &nphg4);
  m_ttree->Branch("nsvtx",                  &nsvtx);

  m_ttree->Branch("G4M_trackid",            &b_G4M_trackid);
  m_ttree->Branch("G4M_nclus",              &b_G4M_nclus);
  m_ttree->Branch("G4M_nclusmvtx",          &b_G4M_nclusmvtx);
  m_ttree->Branch("G4M_nclusintt",          &b_G4M_nclusintt);
  m_ttree->Branch("G4M_nclustpc",           &b_G4M_nclustpc);
  m_ttree->Branch("G4M_nclus_matchrat",     &b_G4M_nclus_matchrat);
  m_ttree->Branch("G4M_nclusmvtx_matchrat", &b_G4M_nclusmvtx_matchrat);
  m_ttree->Branch("G4M_nclusintt_matchrat", &b_G4M_nclusintt_matchrat);
  m_ttree->Branch("G4M_nclustpc_matchrat",  &b_G4M_nclustpc_matchrat);
  m_ttree->Branch("G4M_pt",                 &b_G4M_pt);
  m_ttree->Branch("G4M_phi",                &b_G4M_phi);
  m_ttree->Branch("G4M_eta",                &b_G4M_eta);
  m_ttree->Branch("SvM_trackid",            &b_SvM_trackid);
  m_ttree->Branch("SvM_nclus",              &b_SvM_nclus);
  m_ttree->Branch("SvM_nclusmvtx",          &b_SvM_nclusmvtx);
  m_ttree->Branch("SvM_nclusintt",          &b_SvM_nclusintt);
  m_ttree->Branch("SvM_nclustpc",           &b_SvM_nclustpc);
  m_ttree->Branch("SvM_nclus_matchrat",     &b_SvM_nclus_matchrat);
  m_ttree->Branch("SvM_nclusmvtx_matchrat", &b_SvM_nclusmvtx_matchrat);
  m_ttree->Branch("SvM_nclusintt_matchrat", &b_SvM_nclusintt_matchrat);
  m_ttree->Branch("SvM_nclustpc_matchrat",  &b_SvM_nclustpc_matchrat);
  m_ttree->Branch("SvM_pt",                 &b_SvM_pt);
  m_ttree->Branch("SvM_phi",                &b_SvM_phi);
  m_ttree->Branch("SvM_eta",                &b_SvM_eta);
  if (m_fill_clusters) {
    m_ttree->Branch("clusM_i0",             &b_clusM_i0);
    m_ttree->Branch("clusM_i1",             &b_clusM_i1);
    m_ttree->Branch("clusM_layer",          &b_clusM_layer);
    m_ttree->Branch("clusM_x",              &b_clusM_x);
    m_ttree->Branch("clusM_y",              &b_clusM_y);
    m_ttree->Branch("clusM_z",              &b_clusM_z);
    m_ttree->Branch("clusM_r",              &b_clusM_r);
    m_ttree->Branch("G4M_clusU_i0",         &b_G4M_clusU_i0);
    m_ttree->Branch("G4M_clusU_i1",         &b_G4M_clusU_i1);
    m_ttree->Branch("G4M_clusU_layer",      &b_G4M_clusU_layer);
    m_ttree->Branch("G4M_clusU_x",          &b_G4M_clusU_x);
    m_ttree->Branch("G4M_clusU_y",          &b_G4M_clusU_y);
    m_ttree->Branch("G4M_clusU_z",          &b_G4M_clusU_z);
    m_ttree->Branch("G4M_clusU_r",          &b_G4M_clusU_r);
    m_ttree->Branch("SvM_clusU_i0",         &b_SvM_clusU_i0);
    m_ttree->Branch("SvM_clusU_i1",         &b_SvM_clusU_i1);
    m_ttree->Branch("SvM_clusU_layer",      &b_SvM_clusU_layer);
    m_ttree->Branch("SvM_clusU_x",          &b_SvM_clusU_x);
    m_ttree->Branch("SvM_clusU_y",          &b_SvM_clusU_y);
    m_ttree->Branch("SvM_clusU_z",          &b_SvM_clusU_z);
    m_ttree->Branch("SvM_clusU_r",          &b_SvM_clusU_r);
  }
  m_ttree->Branch("G4U_trackid",            &b_G4U_trackid);
  m_ttree->Branch("G4U_nclus",              &b_G4U_nclus);
  m_ttree->Branch("G4U_nclusmvtx",          &b_G4U_nclusmvtx);
  m_ttree->Branch("G4U_nclusintt",          &b_G4U_nclusintt);
  m_ttree->Branch("G4U_nclustpc",           &b_G4U_nclustpc);
  m_ttree->Branch("G4U_pt",                 &b_G4U_pt);
  m_ttree->Branch("G4U_phi",                &b_G4U_phi);
  m_ttree->Branch("G4U_eta",                &b_G4U_eta);
  if (m_fill_SvU) {
    m_ttree->Branch("SvU_trackid",          &b_SvU_trackid);
    m_ttree->Branch("SvU_nclus",            &b_SvU_nclus);
    m_ttree->Branch("SvU_nclusmvtx",        &b_SvU_nclusmvtx);
    m_ttree->Branch("SvU_nclusintt",        &b_SvU_nclusintt);
    m_ttree->Branch("SvU_nclustpc",         &b_SvU_nclustpc);
    m_ttree->Branch("SvU_pt",               &b_SvU_pt);
    m_ttree->Branch("SvU_phi",              &b_SvU_phi);
    m_ttree->Branch("SvU_eta",              &b_SvU_eta);
  }
  if (m_fill_clusters) {
    m_ttree->Branch("G4U_clusU_i0",         &b_G4U_clusU_i0);
    m_ttree->Branch("G4U_clusU_i1",         &b_G4U_clusU_i1);
    m_ttree->Branch("G4U_clusU_layer",      &b_G4U_clusU_layer);
    m_ttree->Branch("G4U_clusU_x",          &b_G4U_clusU_x);
    m_ttree->Branch("G4U_clusU_y",          &b_G4U_clusU_y);
    m_ttree->Branch("G4U_clusU_z",          &b_G4U_clusU_z);
    m_ttree->Branch("G4U_clusU_r",          &b_G4U_clusU_r);
    if (m_fill_SvU) {
      m_ttree->Branch("SvU_clusU_i0",       &b_SvU_clusU_i0);
      m_ttree->Branch("SvU_clusU_i1",       &b_SvU_clusU_i1);
      m_ttree->Branch("SvU_clusU_layer",    &b_SvU_clusU_layer);
      m_ttree->Branch("SvU_clusU_x",        &b_SvU_clusU_x);
      m_ttree->Branch("SvU_clusU_y",        &b_SvU_clusU_y);
      m_ttree->Branch("SvU_clusU_z",        &b_SvU_clusU_z);
      m_ttree->Branch("SvU_clusU_r",        &b_SvU_clusU_r);
    }
  }


}

//____________________________________________________________________________..
FillTruthRecoMatchTree::~FillTruthRecoMatchTree()
{
}

//____________________________________________________________________________..
int FillTruthRecoMatchTree::Init(PHCompositeNode *topNode)
{
  if (Verbosity()>1) {
    std::cout << " Beginning FillTruthRecoMatchTree " << std::endl;
    topNode->print();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int FillTruthRecoMatchTree::InitRun(PHCompositeNode *topNode)
{
  auto init_status = m_cluster_comp.init(topNode);
  if (init_status == Fun4AllReturnCodes::ABORTRUN) return init_status;

  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int FillTruthRecoMatchTree::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cout << PHWHERE << " DST node is missing, quitting" << std::endl;
    std::cerr << PHWHERE << " DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in FillTruthRecoMatchTree::createNodes");
  }

  m_EmbRecoMatchContainer = findNode::getClass<EmbRecoMatchContainer>(topNode,"TRKR_EMBRECOMATCHCONTAINER");
  if (!m_EmbRecoMatchContainer) {
    std::cout << PHWHERE << " Cannot find node TRKR_EMBRECOMATCHCONTAINER on node tree; quitting " << std::endl;
    std::cerr << PHWHERE << " Cannot find node TRKR_EMBRECOMATCHCONTAINER on node tree; quitting " << std::endl;
    throw std::runtime_error(" Cannot find node TRKR_EMBRECOMATCHCONTAINER on node tree; quitting");
  }

  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));
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
      <<"\"TruthRecoTrackMatching\" module." << std::endl;
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

int FillTruthRecoMatchTree::process_event(PHCompositeNode * /*topNode*/)
{

  if (Verbosity()>5) cout << " FillTruthRecoMatchTree::process_event() " << endl;

  // fill in the event data
  ++nevent; 
  nphg4 = m_TrkrTruthTrackContainer ->getMap().size();
  nsvtx = m_SvtxTrackMap            ->size();
  ntrackmatches = m_EmbRecoMatchContainer->getMatches().size();
  // get centrality later...
  
  // fill in pixel widths on truth tracks
  
  for (auto hitsetkey : m_cluscntr.get_PHG4_clusters()->getHitSetKeys()) {
    float layer = (float) TrkrDefs::getLayer(hitsetkey);
    auto range = m_cluscntr.get_PHG4_clusters()->getClusters(hitsetkey);
    for (auto iter = range.first; iter != range.second; ++iter) {
      auto& cluster = iter->second;
      h2_G4_nPixelsPhi ->Fill( (float)cluster->getPhiSize(), layer );
      h2_G4_nPixelsZ   ->Fill( (float)cluster->getZSize(),   layer );
    }
  }
  // fill in pixel widths on reco tracks
  for (auto hitsetkey : m_cluscntr.get_SVTX_clusters()->getHitSetKeys()) {
    float layer = (float) TrkrDefs::getLayer(hitsetkey);
    auto range = m_cluscntr.get_SVTX_clusters()->getClusters(hitsetkey);
    for (auto iter = range.first; iter != range.second; ++iter) {
      auto& cluster = iter->second;
      h2_Sv_nPixelsPhi ->Fill( (float)cluster->getPhiSize(), layer );
      h2_Sv_nPixelsZ   ->Fill( (float)cluster->getZSize(), layer );
    }
  }

  nphg4_part = 0;
	const auto range = m_PHG4TruthInfoContainer->GetPrimaryParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
      iter != range.second; ++iter)
  { nphg4_part++; }



  // unmatched tracks are only entered once
  // matches can repeat a given svtx or phg4 track, depending on the 
  // parameters in teh matching in filltruthrecomatchtree
  //
  // (1) fill unmatched phg4
  // (2) fill unmatched svtx
  // (3) fill matched phg4 and svtx
  clear_vectors(); 
  int index_G4M_clusU   {0};
  int index_SvM_clusU   {0};
  int i_matched   {0};
  int index_G4U_clusU   {0};
  int index_SvU_clusU   {0};
  
  if (Verbosity() > 2) std::cout << " getting" << (int) m_EmbRecoMatchContainer->getMatches().size() << std::endl;
  for (auto& match : m_EmbRecoMatchContainer->getMatches()) {

    unsigned int g4_trkid = match->idTruthTrack();
    int sv_trkid = match->idRecoTrack();

    auto g4trk = m_TrkrTruthTrackContainer->getTruthTrack(g4_trkid);
    auto svtrk = m_SvtxTrackMap->get(sv_trkid);

    b_G4M_trackid .push_back(g4_trkid);
    b_G4M_pt      .push_back(g4trk->getPt());
    b_G4M_eta     .push_back(g4trk->getPseudoRapidity());
    b_G4M_phi     .push_back(g4trk->getPhi());

    b_SvM_trackid .push_back(sv_trkid);
    b_SvM_pt      .push_back(svtrk->get_pt());
    b_SvM_eta     .push_back(svtrk->get_eta());
    b_SvM_phi     .push_back(svtrk->get_phi());

    m_cluscntr.addClusKeys( g4trk );
    m_cluscntr.addClusKeys( svtrk );
    m_cluscntr.find_matches();

    auto cnt = m_cluscntr.phg4_cntclus();
    b_G4M_nclus          .push_back( cnt[4] );
    b_G4M_nclusmvtx      .push_back( cnt[0] );
    b_G4M_nclusintt      .push_back( cnt[1] );
    b_G4M_nclustpc       .push_back( cnt[2] );
    auto cnt_match = m_cluscntr.phg4_cnt_matchedclus();
    b_G4M_nclus_matchrat     .push_back( (float)cnt_match[4] / cnt[4] );
    b_G4M_nclusmvtx_matchrat .push_back( (float)cnt_match[0] / cnt[0] );
    b_G4M_nclusintt_matchrat .push_back( (float)cnt_match[1] / cnt[1] );
    b_G4M_nclustpc_matchrat  .push_back( (float)cnt_match[2] / cnt[2] );

    cnt = m_cluscntr.svtx_cntclus();
    b_SvM_nclus              .push_back( cnt[4] );
    b_SvM_nclusmvtx          .push_back( cnt[0] );
    b_SvM_nclusintt          .push_back( cnt[1] );
    b_SvM_nclustpc           .push_back( cnt[2] );
    cnt_match = m_cluscntr.svtx_cnt_matchedclus();
    b_SvM_nclus_matchrat     .push_back( (float)cnt_match[4] / cnt[4] );
    b_SvM_nclusmvtx_matchrat .push_back( (float)cnt_match[0] / cnt[0] );
    b_SvM_nclusintt_matchrat .push_back( (float)cnt_match[1] / cnt[1] );
    b_SvM_nclustpc_matchrat  .push_back( (float)cnt_match[2] / cnt[2] );

    if (m_fill_clusters) {
      // clusters only in G4 matched tracks (i.e. by definition, unmatched clusters)
      auto clusters = m_cluscntr.phg4_clusloc_unmatched();
      b_G4M_clusU_i0.push_back(index_G4M_clusU);
      index_G4M_clusU += clusters.size();
      b_G4M_clusU_i1.push_back(index_G4M_clusU);
      for (auto& loc : clusters) {
        b_G4M_clusU_layer.push_back(loc.first);
        b_G4M_clusU_x.push_back(loc.second[0]);
        b_G4M_clusU_y.push_back(loc.second[1]);
        b_G4M_clusU_z.push_back(loc.second[2]);
        b_G4M_clusU_r.push_back(pow(pow(loc.second[0],2.)+pow(loc.second[1],2.),0.5));
      }

      // Svtx only (i.e. unmatched) clusters
      clusters = m_cluscntr.svtx_clusloc_unmatched();
      b_SvM_clusU_i0.push_back(index_SvM_clusU);
      index_SvM_clusU += clusters.size();
      b_SvM_clusU_i1.push_back(index_SvM_clusU);
      for (auto& loc : clusters) {
        b_SvM_clusU_layer.push_back(loc.first);
        b_SvM_clusU_x.push_back(loc.second[0]);
        b_SvM_clusU_y.push_back(loc.second[1]);
        b_SvM_clusU_z.push_back(loc.second[2]);
        b_SvM_clusU_r.push_back(pow(pow(loc.second[0],2.)+pow(loc.second[1],2.),0.5));
      }

      // Matched clusters
      clusters = m_cluscntr.clusloc_matched();
      b_clusM_i0.push_back(i_matched);
      i_matched += clusters.size();
      b_clusM_i1.push_back(i_matched);
      for (auto& loc : clusters) {
        b_clusM_layer.push_back(loc.first);
        b_clusM_x.push_back(loc.second[0]);
        b_clusM_y.push_back(loc.second[1]);
        b_clusM_z.push_back(loc.second[2]);
        b_clusM_r.push_back(pow(pow(loc.second[0],2.)+pow(loc.second[1],2.),0.5));
      }
    }
  }

  // fill in un-matched PHG4 tracks
  for (auto& g4_trkid : m_EmbRecoMatchContainer->ids_TruthUnmatched()) {
    auto g4trk = m_TrkrTruthTrackContainer->getTruthTrack(g4_trkid);
    b_G4U_trackid .push_back(g4_trkid);
    b_G4U_pt      .push_back(g4trk->getPt());
    b_G4U_eta     .push_back(g4trk->getPseudoRapidity());
    b_G4U_phi     .push_back(g4trk->getPhi());

    m_cluscntr.addClusKeys( g4trk );
    auto cnt = m_cluscntr.phg4_cntclus();
    b_G4U_nclus          .push_back( cnt[4] );
    b_G4U_nclusmvtx      .push_back( cnt[0] );
    b_G4U_nclusintt      .push_back( cnt[1] );
    b_G4U_nclustpc       .push_back( cnt[2] );
    if (m_fill_clusters) {
      auto clusters = m_cluscntr.phg4_clusloc_all();
      b_G4U_clusU_i0.push_back(index_G4U_clusU);
      index_G4U_clusU += clusters.size();
      b_G4U_clusU_i1.push_back(index_G4U_clusU);
      for (auto& loc : clusters) {
        b_G4U_clusU_layer.push_back(loc.first);
        b_G4U_clusU_x.push_back(loc.second[0]);
        b_G4U_clusU_y.push_back(loc.second[1]);
        b_G4U_clusU_z.push_back(loc.second[2]);
        b_G4U_clusU_r.push_back(pow(pow(loc.second[0],2.)+pow(loc.second[1],2.),0.5));
      }
    }
  }

  if (m_fill_SvU) {
    for (auto sv_trkid : G4Eval::unmatchedSvtxTrkIds(m_EmbRecoMatchContainer, m_SvtxTrackMap)) {
      auto svtrk = m_SvtxTrackMap->get(sv_trkid);
      b_SvU_trackid .push_back(sv_trkid);
      b_SvU_pt      .push_back(svtrk->get_pt());
      b_SvU_eta     .push_back(svtrk->get_eta());
      b_SvU_phi     .push_back(svtrk->get_phi());
      m_cluscntr.addClusKeys( svtrk );
      auto cnt = m_cluscntr.svtx_cntclus();
      b_SvU_nclus     .push_back( cnt[4] );
      b_SvU_nclusmvtx .push_back( cnt[0] );
      b_SvU_nclusintt .push_back( cnt[1] );
      b_SvU_nclustpc  .push_back( cnt[2] );
    }
    if (m_fill_clusters) {
      auto clusters = m_cluscntr.svtx_clusloc_all();
      b_SvU_clusU_i0.push_back (index_SvU_clusU);
      index_SvU_clusU += clusters.size();
      b_SvU_clusU_i1.push_back (index_SvU_clusU);
      for (auto& loc : clusters) {
        b_SvU_clusU_layer.push_back(loc.first);
        b_SvU_clusU_x.push_back(loc.second[0]);
        b_SvU_clusU_y.push_back(loc.second[1]);
        b_SvU_clusU_z.push_back(loc.second[2]);
        b_SvU_clusU_r.push_back(pow(pow(loc.second[0],2.)+pow(loc.second[1],2.),0.5));
      }
    }
  }

  m_ttree->Fill();
  clear_vectors();
  if (Verbosity()>100) print_mvtx_diagnostics();
  return Fun4AllReturnCodes::EVENT_OK;
}

void FillTruthRecoMatchTree::print_mvtx_diagnostics() {
  std::cout << "To do: "
    << " (1)  number of truth tracks and total number of mvtx and ratio " << std::endl
    << " (2)  ditto for reco tracks " << std::endl;

  double n_PHG4_tracks = m_TrkrTruthTrackContainer->getMap().size();
  // count how many mvtx clusters in the phg4 tracks
  double n_in_PHG4_tracks {0.};
  for (auto& pair : m_TrkrTruthTrackContainer->getMap()) {
    m_cluscntr.addClusKeys( pair.second );
    n_in_PHG4_tracks += m_cluscntr.phg4_cntclus()[0];
  }
  // count how mant mvtx clusters in truth container (should be the same)
  double n_in_PHG4_clusters {0.};
  for (auto hitsetkey : m_cluscntr.get_PHG4_clusters()->getHitSetKeys()) {
    if (TrkrDefs::getLayer(hitsetkey) > 2) continue;
    auto range = m_cluscntr.get_PHG4_clusters()->getClusters(hitsetkey);
    for (auto r = range.first; r != range.second; ++r)  {
      n_in_PHG4_clusters += 1.;
    }
  }

  // count how many svtx tracks
  double n_SVTX_tracks = m_SvtxTrackMap->size();
  // count how many mvtx clusters in svtx tracks
  double n_in_SVTX_tracks {0.};
  for (auto entry = m_SvtxTrackMap->begin(); entry != m_SvtxTrackMap->end(); ++entry) {
    m_cluscntr.addClusKeys(entry->second);
    n_in_SVTX_tracks += m_cluscntr.svtx_cntclus()[0];
  }
  // count how many mvtx are total in the container
  double n_in_SVTX_clusters {0.};
  for (auto hitsetkey : m_cluscntr.get_SVTX_clusters()->getHitSetKeys()) {
    if (TrkrDefs::getLayer(hitsetkey) > 2) continue;
    auto range = m_cluscntr.get_SVTX_clusters()->getClusters(hitsetkey);
    for (auto r = range.first; r != range.second; ++r)  {
      n_in_SVTX_clusters += 1.;
    }
  }

  std::cout << Form("MVTX"
      "\nPHG4:  Tracks(%.0f)   Clusters In tracks(%.0f)   Total (%.0f)"
      "\n       ave. per track: %6.3f   ratio in all tracks: %6.2f"
      , n_PHG4_tracks, n_in_PHG4_tracks, n_in_PHG4_clusters
      , (n_in_PHG4_tracks/n_PHG4_tracks), (n_in_PHG4_tracks/n_in_PHG4_clusters))
    << std::endl;
  std::cout << Form(
      "\nSVTX:  Tracks(%.0f)   Clusters In tracks(%.0f)   Total (%.0f)"
      "\n       ave. per track: %6.3f   ratio in all tracks: %6.2f"
      , n_SVTX_tracks, n_in_SVTX_tracks, n_in_SVTX_clusters
      , (n_in_SVTX_tracks/n_SVTX_tracks), (n_in_SVTX_tracks/n_in_SVTX_clusters))
    << std::endl;

}

int FillTruthRecoMatchTree::End(PHCompositeNode *)
{
  if (Verbosity()>2) std::cout << PHWHERE << ": ending FillTruthRecoMatchTree" << std::endl;
	PHTFileServer::get().cd(m_outfile_name);

  h2_G4_nPixelsPhi ->Write();
  h2_G4_nPixelsZ   ->Write();
  h2_Sv_nPixelsPhi ->Write();
  h2_Sv_nPixelsZ   ->Write();

  m_ttree->Write();
  return Fun4AllReturnCodes::EVENT_OK;
}

void FillTruthRecoMatchTree::clear_vectors() {
    // Tracks and clustes

    b_G4M_trackid            .clear();
    b_G4M_nclus              .clear();
    b_G4M_nclusmvtx          .clear();
    b_G4M_nclusintt          .clear();
    b_G4M_nclustpc           .clear();
    b_G4M_nclus_matchrat     .clear();
    b_G4M_nclusmvtx_matchrat .clear();
    b_G4M_nclusintt_matchrat .clear();
    b_G4M_nclustpc_matchrat  .clear();
    b_G4M_pt                 .clear();
    b_G4M_phi                .clear();
    b_G4M_eta                .clear();
    b_SvM_trackid            .clear();
    b_SvM_nclus              .clear();
    b_SvM_nclusmvtx          .clear();
    b_SvM_nclusintt          .clear();
    b_SvM_nclustpc           .clear();
    b_SvM_nclus_matchrat     .clear();
    b_SvM_nclusmvtx_matchrat .clear();
    b_SvM_nclusintt_matchrat .clear();
    b_SvM_nclustpc_matchrat  .clear();
    b_SvM_pt                 .clear();
    b_SvM_phi                .clear();
    b_SvM_eta                .clear();
    if (m_fill_clusters) {
      b_clusM_i0               .clear();
      b_clusM_i1               .clear();
      b_clusM_layer            .clear();
      b_clusM_x                .clear();
      b_clusM_y                .clear();
      b_clusM_z                .clear();
      b_clusM_r                .clear();
      b_G4M_clusU_i0           .clear();
      b_G4M_clusU_i1           .clear();
      b_G4M_clusU_layer        .clear();
      b_G4M_clusU_x            .clear();
      b_G4M_clusU_y            .clear();
      b_G4M_clusU_z            .clear();
      b_G4M_clusU_r            .clear();
      b_SvM_clusU_i0           .clear();
      b_SvM_clusU_i1           .clear();
      b_SvM_clusU_layer        .clear();
      b_SvM_clusU_x            .clear();
      b_SvM_clusU_y            .clear();
      b_SvM_clusU_z            .clear();
      b_SvM_clusU_r            .clear();
    }
    b_G4U_trackid            .clear();
    b_G4U_nclus              .clear();
    b_G4U_nclusmvtx          .clear();
    b_G4U_nclusintt          .clear();
    b_G4U_nclustpc           .clear();
    b_G4U_pt                 .clear();
    b_G4U_phi                .clear();
    b_G4U_eta                .clear();
    if (m_fill_SvU) {
      b_SvU_trackid            .clear();
      b_SvU_nclus              .clear();
      b_SvU_nclusmvtx          .clear();
      b_SvU_nclusintt          .clear();
      b_SvU_nclustpc           .clear();
      b_SvU_pt                 .clear();
      b_SvU_phi                .clear();
      b_SvU_eta                .clear();
    }
    if (m_fill_clusters) {
      b_G4U_clusU_i0           .clear();
      b_G4U_clusU_i1           .clear();
      b_G4U_clusU_layer        .clear();
      b_G4U_clusU_x            .clear();
      b_G4U_clusU_y            .clear();
      b_G4U_clusU_z            .clear();
      b_G4U_clusU_r            .clear();
      if (m_fill_SvU) {
        b_SvU_clusU_i0           .clear();
        b_SvU_clusU_i1           .clear();
        b_SvU_clusU_layer        .clear();
        b_SvU_clusU_x            .clear();
        b_SvU_clusU_y            .clear();
        b_SvU_clusU_z            .clear();
        b_SvU_clusU_r            .clear();
      }
    }



}
