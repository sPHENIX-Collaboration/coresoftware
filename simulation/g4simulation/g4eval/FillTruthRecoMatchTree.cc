#include "FillTruthRecoMatchTree.h"

#include <TFile.h>
#include <TTree.h>
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
#include <trackbase_historic/SvtxPHG4ParticleMap_v1.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

using std::cout;
using std::endl;
//____________________________________________________________________________..
FillTruthRecoMatchTree::FillTruthRecoMatchTree(const std::string &name
    , bool _fill_clusters
    , const std::string _tfile_name)
  : SubsysReco(name), m_fill_clusters { _fill_clusters }
{
  m_tfile                                   = new TFile(_tfile_name.c_str(), "recreate");
  m_tfile->cd();
  m_ttree->Branch("event",                  &nevent);
  m_ttree->Branch("nphg4_part",             &nphg4_part);
  m_ttree->Branch("centrality",             &centrality);
  m_ttree->Branch("ntrackmatches",          &nmatchtracks);
  m_ttree->Branch("nphg4",                  &nphg4);
  m_ttree->Branch("nsvtx",                  &nsvtx);

  m_ttree->Branch("g4M_trackid",            &b_g4M_trackid);
  m_ttree->Branch("g4M_nclus",              &b_g4M_nclus);
  m_ttree->Branch("g4M_nclusmvtx",          &b_g4M_nclusmvtx);
  m_ttree->Branch("g4M_nclusintt",          &b_g4M_nclusintt);
  m_ttree->Branch("g4M_nclustpc",           &b_g4M_nclustpc);
  m_ttree->Branch("g4M_nclusmvtx_matchrat", &b_g4M_nclusmvtx_matchrat);
  m_ttree->Branch("g4M_nclusintt_matchrat", &b_g4M_nclusintt_matchrat);
  m_ttree->Branch("g4M_nclustpc_matchrat",  &b_g4M_nclustpc_matchrat);
  m_ttree->Branch("g4M_pt",                 &b_g4M_pt);
  m_ttree->Branch("g4M_px",                 &b_g4M_px);
  m_ttree->Branch("g4M_py",                 &b_g4M_py);
  m_ttree->Branch("g4M_pz",                 &b_g4M_pz);
  m_ttree->Branch("g4M_clusM_i0",           &b_g4M_clusM_i0);
  m_ttree->Branch("g4M_clusM_i1",           &b_g4M_clusM_i1);
  m_ttree->Branch("g4M_clusM_layer",        &b_g4M_clusM_layer);
  m_ttree->Branch("g4M_clusM_x",            &b_g4M_clusM_x);
  m_ttree->Branch("g4M_clusM_y",            &b_g4M_clusM_y);
  m_ttree->Branch("g4M_clusM_z",            &b_g4M_clusM_z);
  m_ttree->Branch("g4M_clusM_r",            &b_g4M_clusM_r);
  m_ttree->Branch("g4M_clusU_i0",           &b_g4M_clusU_i0);
  m_ttree->Branch("g4M_clusU_i1",           &b_g4M_clusU_i1);
  m_ttree->Branch("g4M_clusU_layer",        &b_g4M_clusU_layer);
  m_ttree->Branch("g4M_clusU_x",            &b_g4M_clusU_x);
  m_ttree->Branch("g4M_clusU_y",            &b_g4M_clusU_y);
  m_ttree->Branch("g4M_clusU_z",            &b_g4M_clusU_z);
  m_ttree->Branch("g4M_clusU_r",            &b_g4M_clusU_r);
  m_ttree->Branch("svM_trackid",            &b_svM_trackid);
  m_ttree->Branch("svM_nclus",              &b_svM_nclus);
  m_ttree->Branch("svM_nclusmvtx",          &b_svM_nclusmvtx);
  m_ttree->Branch("svM_nclusintt",          &b_svM_nclusintt);
  m_ttree->Branch("svM_nclustpc",           &b_svM_nclustpc);
  m_ttree->Branch("svM_nclusmvtx_matchrat", &b_svM_nclusmvtx_matchrat);
  m_ttree->Branch("svM_nclusintt_matchrat", &b_svM_nclusintt_matchrat);
  m_ttree->Branch("svM_nclustpc_matchrat",  &b_svM_nclustpc_matchrat);
  m_ttree->Branch("svM_pt",                 &b_svM_pt);
  m_ttree->Branch("svM_px",                 &b_svM_px);
  m_ttree->Branch("svM_py",                 &b_svM_py);
  m_ttree->Branch("svM_pz",                 &b_svM_pz);
  m_ttree->Branch("svM_clusM_i0",           &b_svM_clusM_i0);
  m_ttree->Branch("svM_clusM_i1",           &b_svM_clusM_i1);
  m_ttree->Branch("svM_clusM_layer",        &b_svM_clusM_layer);
  m_ttree->Branch("svM_clusM_x",            &b_svM_clusM_x);
  m_ttree->Branch("svM_clusM_y",            &b_svM_clusM_y);
  m_ttree->Branch("svM_clusM_z",            &b_svM_clusM_z);
  m_ttree->Branch("svM_clusM_r",            &b_svM_clusM_r);
  m_ttree->Branch("svM_clusU_i0",           &b_svM_clusU_i0);
  m_ttree->Branch("svM_clusU_i1",           &b_svM_clusU_i1);
  m_ttree->Branch("svM_clusU_layer",        &b_svM_clusU_layer);
  m_ttree->Branch("svM_clusU_x",            &b_svM_clusU_x);
  m_ttree->Branch("svM_clusU_y",            &b_svM_clusU_y);
  m_ttree->Branch("svM_clusU_z",            &b_svM_clusU_z);
  m_ttree->Branch("svM_clusU_r",            &b_svM_clusU_r);
  m_ttree->Branch("nclus_match",            &b_nclus_match);
  m_ttree->Branch("g4svmatch_nclusmvtx",    &b_g4svmatch_nclusmvtx);
  m_ttree->Branch("g4svmatch_nclusintt",    &b_g4svmatch_nclusintt);
  m_ttree->Branch("g4svmatch_nclustpc",     &b_g4svmatch_nclustpc);
  m_ttree->Branch("g4U_trackid",            &b_g4U_trackid);
  m_ttree->Branch("g4U_nclus",              &b_g4U_nclus);
  m_ttree->Branch("g4U_nclusmvtx",          &b_g4U_nclusmvtx);
  m_ttree->Branch("g4U_nclusintt",          &b_g4U_nclusintt);
  m_ttree->Branch("g4U_nclustpc",           &b_g4U_nclustpc);
  m_ttree->Branch("g4U_pt",                 &b_g4U_pt);
  m_ttree->Branch("g4U_px",                 &b_g4U_px);
  m_ttree->Branch("g4U_py",                 &b_g4U_py);
  m_ttree->Branch("g4U_pz",                 &b_g4U_pz);
  m_ttree->Branch("g4U_clusU_i0",           &b_g4U_clusU_i0);
  m_ttree->Branch("g4U_clusU_i1",           &b_g4U_clusU_i1);
  m_ttree->Branch("g4U_clusU_layer",        &b_g4U_clusU_layer);
  m_ttree->Branch("g4U_clusU_x",            &b_g4U_clusU_x);
  m_ttree->Branch("g4U_clusU_y",            &b_g4U_clusU_y);
  m_ttree->Branch("g4U_clusU_z",            &b_g4U_clusU_z);
  m_ttree->Branch("g4U_clusU_r",            &b_g4U_clusU_r);
  m_ttree->Branch("svU_trackid",            &b_svU_trackid);
  m_ttree->Branch("svU_nclus",              &b_svU_nclus);
  m_ttree->Branch("svU_nclusmvtx",          &b_svU_nclusmvtx);
  m_ttree->Branch("svU_nclusintt",          &b_svU_nclusintt);
  m_ttree->Branch("svU_nclustpc",           &b_svU_nclustpc);
  m_ttree->Branch("svU_pt",                 &b_svU_pt);
  m_ttree->Branch("svU_px",                 &b_svU_px);
  m_ttree->Branch("svU_py",                 &b_svU_py);
  m_ttree->Branch("svU_pz",                 &b_svU_pz);
  m_ttree->Branch("svU_clusU_i0",           &b_svU_clusU_i0);
  m_ttree->Branch("svU_clusU_i1",           &b_svU_clusU_i1);
  m_ttree->Branch("svU_clusU_layer",        &b_svU_clusU_layer);
  m_ttree->Branch("svU_clusU_x",            &b_svU_clusU_x);
  m_ttree->Branch("svU_clusU_y",            &b_svU_clusU_y);
  m_ttree->Branch("svU_clusU_z",            &b_svU_clusU_z);
  m_ttree->Branch("svU_clusU_r",            &b_svU_clusU_r);

}

//____________________________________________________________________________..
FillTruthRecoMatchTree::~FillTruthRecoMatchTree()
{
}

//____________________________________________________________________________..
int FillTruthRecoMatchTree::Init(PHCompositeNode *topNode)
{
  if (Verbosity()>1) topNode->print();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int FillTruthRecoMatchTree::InitRun(PHCompositeNode *topNode)
{
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


  m_TruthClusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, 
      "TRKR_TRUTHCLUSTERCONTAINER");
  if (!m_TruthClusterContainer)
  {
    std::cout << "Could not locate TRKR_TRUTHCLUSTERCONTAINER node when running "
      << "\"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_RecoClusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_RecoClusterContainer)
  {
    std::cout << "Could not locate TRKR_CLUSTER node when running "
      << "\"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

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
  for (auto& match : m_EmbRecoMatchContainer->getMatches()) {
    if (match) cout << "waste" << endl;
    // every match is a unique match from truth to embedded particles
    /* const int            gtrackID = match->idTruthTrack()   ; */
    /* const unsigned int   id_reco  = match->idRecoTrack()    ; */
    /* const unsigned short n_match  = match->nClustersMatched() ; */
    /* const unsigned short n_truth  = match->nClustersTruth() ; */
    /* const unsigned short n_reco   = match->nClustersReco()  ; */

    /* // fill the map_TtoR */
    /* /1* PHG4ParticleSvtxMap::WeightedRecoTrackMap*  entry_TtoR; *1/ */
    /* if (map_TtoR.find(gtrackID) == map_TtoR.end()) { */
    /*   map_TtoR[gtrackID] = PHG4ParticleSvtxMap::WeightedRecoTrackMap{}; */
    /* } */
    /* auto& entry_TtoR = map_TtoR[gtrackID]; */
    /* float weight_TtoR = (float)n_match + (float)n_truth/100.; */
    /* if (entry_TtoR.find(weight_TtoR) == entry_TtoR.end()) { */
    /*   entry_TtoR[weight_TtoR] = { id_reco }; // i.e. std::set<unsigned int> { id_reco }; */
    /* } else { */
    /*   entry_TtoR[weight_TtoR].insert(id_reco); */
    /* } */

    /* // fill the map_RtoT */
    /* /1* SvtxPHG4ParticleMap::WeightedTruthTrackMap*  entry_RtoT; *1/ */
    /* if (map_RtoT.find(id_reco) == map_RtoT.end()) { */
    /*   map_RtoT[id_reco] = SvtxPHG4ParticleMap::WeightedTruthTrackMap{}; */
    /* } */ 
    /* auto& entry_RtoT = map_RtoT[id_reco]; */
    /* float weight_RtoT = (float)n_match + (float)n_reco/100.; */
    /* if (entry_RtoT.find(weight_RtoT) == entry_RtoT.end()) { */
    /*     entry_RtoT[weight_RtoT] = { gtrackID }; // i.e. std::set<int> { gtrackID } */
    /* } else { */
    /*   entry_RtoT[weight_RtoT].insert(gtrackID); */
    /* } */

    /* if (Verbosity() > 20) { */
    /*   printf("EmbRecoMatch: gtrackID(%2i) id_reco(%2i) nclusters:match(%i),gtrack(%2i),reco(%2i)\n", */
    /*       gtrackID, (int)id_reco, (int)n_match, (int)n_truth, (int)n_reco); */
    /*   printf("   -> in SvtxPHG4ParticleMap {id_reco->{weight->id_true}} = {%2i->{%5.2f->%2i}}\n", */ 
    /*       (int)id_reco,  weight_RtoT, (int)gtrackID); */
    /*   printf("   -> in PHG4ParticleSvtxMap {id_true->{weight->id_reco}} = {%2i->{%5.2f->%2i}}\n", */ 
    /*      gtrackID, weight_TtoR, (int)id_reco ); */
  }

  m_ttree->Fill();
  clear_vectors();
  return Fun4AllReturnCodes::EVENT_OK;
}

int FillTruthRecoMatchTree::End(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void FillTruthRecoMatchTree::clear_vectors() {
    // Tracks and clustes
    //
    // lables:
    //  tracks:
    //    g4M : phg4track matched
    //    svM : svtx_track matched
    //    g4U : phg4track not-matched
    //    svU : svtx_track not-matched
    //  clusters:
    //    clusM : matched
    //    clusU : unmatched
    b_g4M_trackid            .clear();
    b_g4M_nclus              .clear();
    b_g4M_nclusmvtx          .clear();
    b_g4M_nclusintt          .clear();
    b_g4M_nclustpc           .clear();
    b_g4M_nclusmvtx_matchrat .clear();
    b_g4M_nclusintt_matchrat .clear();
    b_g4M_nclustpc_matchrat  .clear();
    b_g4M_pt                 .clear();
    b_g4M_px                 .clear();
    b_g4M_py                 .clear();
    b_g4M_pz                 .clear();
    b_g4M_clusM_i0           .clear();
    b_g4M_clusM_i1           .clear();
    b_g4M_clusM_layer        .clear();
    b_g4M_clusM_x            .clear();
    b_g4M_clusM_y            .clear();
    b_g4M_clusM_z            .clear();
    b_g4M_clusM_r            .clear();
    b_g4M_clusU_i0           .clear();
    b_g4M_clusU_i1           .clear();
    b_g4M_clusU_layer        .clear();
    b_g4M_clusU_x            .clear();
    b_g4M_clusU_y            .clear();
    b_g4M_clusU_z            .clear();
    b_g4M_clusU_r            .clear();
    b_svM_trackid            .clear();
    b_svM_nclus              .clear();
    b_svM_nclusmvtx          .clear();
    b_svM_nclusintt          .clear();
    b_svM_nclustpc           .clear();
    b_svM_nclusmvtx_matchrat .clear();
    b_svM_nclusintt_matchrat .clear();
    b_svM_nclustpc_matchrat  .clear();
    b_svM_pt                 .clear();
    b_svM_px                 .clear();
    b_svM_py                 .clear();
    b_svM_pz                 .clear();
    b_svM_clusM_i0           .clear();
    b_svM_clusM_i1           .clear();
    b_svM_clusM_layer        .clear();
    b_svM_clusM_x            .clear();
    b_svM_clusM_y            .clear();
    b_svM_clusM_z            .clear();
    b_svM_clusM_r            .clear();
    b_svM_clusU_i0           .clear();
    b_svM_clusU_i1           .clear();
    b_svM_clusU_layer        .clear();
    b_svM_clusU_x            .clear();
    b_svM_clusU_y            .clear();
    b_svM_clusU_z            .clear();
    b_svM_clusU_r            .clear();
    b_nclus_match            .clear();
    b_g4svmatch_nclusmvtx    .clear();
    b_g4svmatch_nclusintt    .clear();
    b_g4svmatch_nclustpc     .clear();
    b_g4U_trackid            .clear();
    b_g4U_nclus              .clear();
    b_g4U_nclusmvtx          .clear();
    b_g4U_nclusintt          .clear();
    b_g4U_nclustpc           .clear();
    b_g4U_pt                 .clear();
    b_g4U_px                 .clear();
    b_g4U_py                 .clear();
    b_g4U_pz                 .clear();
    b_g4U_clusU_i0           .clear();
    b_g4U_clusU_i1           .clear();
    b_g4U_clusU_layer        .clear();
    b_g4U_clusU_x            .clear();
    b_g4U_clusU_y            .clear();
    b_g4U_clusU_z            .clear();
    b_g4U_clusU_r            .clear();
    b_svU_trackid            .clear();
    b_svU_nclus              .clear();
    b_svU_nclusmvtx          .clear();
    b_svU_nclusintt          .clear();
    b_svU_nclustpc           .clear();
    b_svU_pt                 .clear();
    b_svU_px                 .clear();
    b_svU_py                 .clear();
    b_svU_pz                 .clear();
    b_svU_clusU_i0           .clear();
    b_svU_clusU_i1           .clear();
    b_svU_clusU_layer        .clear();
    b_svU_clusU_x            .clear();
    b_svU_clusU_y            .clear();
    b_svU_clusU_z            .clear();
    b_svU_clusU_r            .clear();
    
}
