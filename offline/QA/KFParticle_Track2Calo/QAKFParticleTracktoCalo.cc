#include "QAKFParticleTracktoCalo.h"

#include <qautils/QAHistManagerDef.h>  // for getHistoManager

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TH2.h>

#include <TSystem.h>

//____________________________________________________________________________..
QAKFParticleTracktoCalo::QAKFParticleTracktoCalo(const std::string& name, const std::string node_name)
  : SubsysReco(name), m_KFParticleNodeName(node_name)
{
}

//____________________________________________________________________________..
int QAKFParticleTracktoCalo::InitRun(PHCompositeNode* topNode)
{  
  
  if(m_verbosity > 10){
    std::cout << " QAKFParticleTracktoCalo: Init run  " << std::endl; 

  }

  // Load geometry node
  if (!EMCalGeo)
  {
    EMCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    if (!EMCalGeo){
      std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "TOWERGEOM_CEMC" << std::endl;
      // return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  
  // Fun4All server 
  auto *se = Fun4AllServer::instance();
  m_runnumber = se->RunNumber();

  // Set up histograms
  SetUpHistos();

  // Return
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void QAKFParticleTracktoCalo::LoadNodes(PHCompositeNode *topNode){

  if(m_verbosity > 10){
    std::cout << " QAKFParticleTracktoCalo: Loading nodes " << std::endl; 
  }

  // Load Track Map
  if (!m_trackMap){ 
    m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
    if (!m_trackMap){
      std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find SvtxTrackMap node " << m_trackMapName << std::endl;      
      //return Fun4AllReturnCodes::ABORTRUN;
    }
  } 
  std::string geoName = "ActsGeometry";
  geometry = findNode::getClass<ActsGeometry>(topNode, geoName);
  if (!geometry)
  { 
    std::cout << "KFParticle detector info: " << geoName << " does not exist" << std::endl;
  } 
  m_dst_clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_dst_clustermap)
  { 
    m_dst_clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_SEED");
    if (!m_dst_clustermap)
    {
      std::cout << "KFParticle detector info: TRKR_CLUSTER does not exist" << std::endl;
    }
  } 

  // Load Clusters
  if (!clustersEM)
  {  
    clustersEM = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");
    if (!clustersEM)
    {
      std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "CLUSTERINFO_CEMC" << std::endl;
    }
  } 
  // Load towers 
  if(!_towersEM)
  { 
    _towersEM = findNode::getClass<TowerInfoContainerv4>(topNode, "TOWERINFO_CALIB_CEMC");
    if(!_towersEM)
    {
      std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "TOWER_CALIB_CEMC" << std::endl;
      //return Fun4AllReturnCodes::ABORTEVENT;
    }
  } 
  //Get KFP container
  std::string kfpContainerNodeName = m_KFParticleNodeName + "_KFParticle_Container";
  if (!m_kfpContainer){ 
    m_kfpContainer = findNode::getClass<KFParticle_Container>(topNode, kfpContainerNodeName);
  } 
  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!m_vertexMap)
  { 
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error : missing SvtxVertexMap" << std::endl;
    // return Fun4AllReturnCodes::ABORTRUN;
  } 
}

//____________________________________________________________________________..
int QAKFParticleTracktoCalo::process_event(PHCompositeNode* topNode)
{
  if(m_verbosity > 10) std::cout << " QAKFParticleTracktoCalo: processing event " << std::endl;
  
  // Load nodes
  LoadNodes(topNode);
  if (!m_kfpContainer){ 
    std::cout << "KFP CONTAINER - Fatal Error - unable to find DST node KFParticle_QA" << std::endl;
    emptyKFPcontainer++;
    return -1;
  }
  else filledKFPcontainer++;
  if(m_verbosity > 10) std::cout << "nParticles in KFP Container:" << m_kfpContainer->size() << std::endl;

  //Loop over the mothers
  for (auto &iter : *m_kfpContainer)
  { std::cout << __LINE__ << std::endl;
    // Setup    // std::vector<std::pair<SvtxTrack*,KFParticle*>> v_daughters={}; // Holds all daughters belonging to a (single) candidate mother
    std::vector<QAKFParticleTracktoCalo::MatchedCluster*> v_matched_clusterstracks = {}; // Holds only duaghters with matches
    SvtxVertex *m_vertex; 

    //Get mother
    KFParticle *mother = iter.second;
    if (!mother){
      if (m_verbosity > 0) std::cout << __PRETTY_FUNCTION__ << " Error: null mother pointer" << std::endl;
      continue;
    } std::cout << __LINE__ << std::endl;

    // Get daughter IDs corresponding to this particle
    std::vector<int> track_ids = mother->DaughterIds();  //KFP container, get mother, get the daughter IDs, should be the same as the track ID        
    
    // Skip daughters -- track_ids.size()=0 for daughters
    if (track_ids.size() < 2) {std::cout << "Bad track_ids.size() < 2..." << std::endl; continue;}
    
    // Get daughters from KFP container
    std::vector<std::pair<SvtxTrack*,KFParticle*>> v_daughters = {};
    v_daughters = GetVectorOfDaughters(m_kfpContainer, m_trackMap, track_ids);

    // Make sure we have the correct number of daughters !! 
    if((int) v_daughters.size() != (int) mother->NDaughters() || (int) v_daughters.size() == 0){
      failedcut_ndaugh++;
      continue;
    }
    
    // Cut on angles between daughters -- CAN ONLY BE APPLIED IF THE DECAY HAS TWO DAUGHTERS!!!!!
    std::pair<float,float> anglesbetweendaughters = {std::numeric_limits<float>::quiet_NaN(),std::numeric_limits<float>::quiet_NaN()};
    if(m_cut_on_angle_between_tracks){
      anglesbetweendaughters = GetEtaAndPhiBetweenDaughters(v_daughters);
      if(!GoodAngleBetweenDaughters(anglesbetweendaughters)){
        failedcut_dAngle++; 
        std::cout << "Bad dAngle: " << anglesbetweendaughters.first << " " << anglesbetweendaughters.second << std::endl; 
        continue;
      }
    }
    
    // Testing -- make sure the vertices are the same !!   
    bool mismatched_vertex=false;
    unsigned int first_id = v_daughters[0].first->get_vertex_id();
    for(const auto& [trk,daughter] : v_daughters){
      if(trk->get_vertex_id() != first_id){
        std::cout << "Event " << m_event << " has candidate with mismatched daughter vertex IDs." << std::endl;
        mismatched_vertex=true;
        break;
      }
    }
    // Printouts for pileup questions
    if(m_verbosity>2 && mismatched_vertex){
      // std::vector<int> v_daughterIDs=mother->DaughterIds();
      // std::cout << "Skipped daughters ";
      // for(const auto& ID : v_daughterIDs){std::cout << ID << " ";}
      // std::cout << "with vertex IDs: ";
      // for(const auto& daughter : v_daughters){std::cout << daughter.first->get_vertex_id() << " ";}
      failedcut_vertex++;
      continue;
    }

    // Get vertex
    m_vertex = m_vertexMap->get(v_daughters[0].first->get_vertex_id());

    // Get matched clusters
    for(const auto& [trk,daughter] : v_daughters){
      // Make sure there is a track
      if(!trk) continue;
      // Get nStates in trackers
      unsigned int nStatesMVTX = 0;
      unsigned int nStatesINTT = 0;
      unsigned int nStatesTPC  = 0;
      GetDetectorHits(trk, nStatesMVTX, nStatesINTT, nStatesTPC);
      // Cut on MVTX hits (optional)
      if(m_cut_on_MaxHitsMVTX && nStatesMVTX > m_MVTX_maxhits) {failedcut_nMVTX++; std::cout << "Bad nMVTX:  " << nStatesMVTX << std::endl; continue;}
      // Get track state
      SvtxTrackState *thisState = nullptr;
      thisState = trk->get_state(caloRadiusEMCal);      
      QAKFParticleTracktoCalo::MatchedCluster *matchedcluster = GetMatchedCluster(thisState, m_vertex);
      if(!matchedcluster->matched) {
        failedcut_nomtch++;
        continue;
      }
      // Get misc. info if the match was found
      matchedcluster->track_charge = trk->get_charge();
      matchedcluster->track_pt_svtx= trk->get_pt();
      matchedcluster->track_p_svtx = trk->get_p();
      matchedcluster->track_eta_svtx= trk->get_eta();
      matchedcluster->track_phi_svtx= trk->get_phi();
      matchedcluster->nStatesINTT= nStatesINTT;
      matchedcluster->nStatesMVTX= nStatesMVTX;
      matchedcluster->nStatesTPC = nStatesTPC;
      v_matched_clusterstracks.push_back(matchedcluster);
    }

    // Fill histograms
    if(v_matched_clusterstracks.size()>0){
      // FillHistos(v_matched_clusterstracks, v_daughters, mother->GetMass(), anglesbetweendaughters);
      FillHistos(v_matched_clusterstracks, mother->GetMass(), anglesbetweendaughters);
    }

  } // End of loop over KFP container for event i

  // End event and return
  EndEvent(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int QAKFParticleTracktoCalo::End(PHCompositeNode * /*topNode*/){ 

  if(m_verbosity > 10){
    std::cout << " QAKFParticleTracktoCalo: Ending " << std::endl;
  }

  // Create out file if not using the histo manager
  if(!use_qa_histomanager){

    // Create Outfile
    std::string m_outfile_namepath = GetOutputFileNameAndPath(m_outfilename);
    std::cout << m_outfilename << std::endl;m_outfile = new TFile(m_outfile_namepath.c_str(), "RECREATE");  // m_outfile = new TFile("TEST.root", "RECREATE");
    
    // Write histos to the .root file
    WriteHistos();
    m_outfile->Write();
    
    // Close outfile  
    m_outfile->Close();

  }
  if(m_verbosity > 10){
    std::cout<< "___________________________________________________________________________" << std::endl;
    std::cout<< "# events with empty KFP containers : " << emptyKFPcontainer  << std::endl;
    std::cout<< "# events with filled KFP containers: " << filledKFPcontainer << std::endl;
    std::cout<< "___________________________________________________________________________" << std::endl;
    std::cout<< "# daughters with failing track cuts..." << std::endl;
    std::cout<< "failedcut_nMVTX:  " << failedcut_nMVTX << std::endl;
    std::cout<< "failedcut_dAngle: " << failedcut_dAngle << std::endl;
    std::cout<< "failedcut_vertex: " << failedcut_vertex << std::endl;
    std::cout<< "failedcut_ndaugh: " << failedcut_ndaugh << std::endl;
    std::cout<< "failedcut_nomtch: " << failedcut_nomtch << std::endl;
    std::cout<< "___________________________________________________________________________" << std::endl;
    std::cout<< "# daughters with failing EMC matching cuts..." << std::endl;
    std::cout<< "failedcut_energy: " << failedcut_energy << std::endl;
    std::cout<< "failedcut_dPhi:   " << failedcut_dPhi << std::endl;
    std::cout<< "failedcut_dZ:     " << failedcut_dZ << std::endl;
    std::cout<< "___________________________________________________________________________" << std::endl;
    
    std::cout << "Number of daughter candidates matched to EMCal clusters: " << number_of_matches << std::endl;
    std::cout<< "___________________________________________________________________________" << std::endl;
  }
  return 0;
}

//____________________________________________________________________________..
void QAKFParticleTracktoCalo::EndEvent(PHCompositeNode*){
  // Get ready for next event
  m_event++;
  clustersEM=nullptr;
  _towersEM=nullptr;
}


//--------------------------------------------------------------------------------
//---------------------------HISTOGRAM FUNCTIONS----------------------------------
//--------------------------------------------------------------------------------
void QAKFParticleTracktoCalo::SetUpHistos(){
  if(m_verbosity > 10){
    std::cout << " QAKFParticleTracktoCalo: Setting up histos " << std::endl;
  }
  
  // QA histo manager
  Fun4AllHistoManager *hm=nullptr;
  if(use_qa_histomanager) {
    hm = QAHistManagerDef::getHistoManager();
    assert(hm);
  }

  // pT
  h_pt_m = new TH1F(Form("%spt_m", histo_prefix.c_str()), "Negative daughter p_{T};p_{T} (GeV/c);Counts", 50, 0, 10);
  h_pt_p = new TH1F(Form("%spt_p", histo_prefix.c_str()), "Positive daughter p_{T};p_{T} (GeV/c);Counts", 50, 0, 10);
  h_pt_m_vs_p = new TH2F(Form("%spt_m_vs_p", histo_prefix.c_str()), "Daughter p_{T} (Both EMCal Cluster Matched); + daughter p_{T} (GeV/c);- daughter p_{T} (GeV/c)", 50, 0, 10, 50, 0, 10);
  // Eta and Phi
  h_eta_vs_phi_m = new TH2F(Form("%seta_vs_phi_m", histo_prefix.c_str()), "Negative daughter #eta vs #varphi;#phi;#eta", 64, -3.2, 3.2, 60, -1.5, 1.5);
  h_eta_vs_phi_p = new TH2F(Form("%seta_vs_phi_p", histo_prefix.c_str()), "Positive daughter #eta vs #varphi;#phi;#eta", 64, -3.2, 3.2, 60, -1.5, 1.5);
  // E/p
  h_eop_m = new TH1F(Form("%seop_m", histo_prefix.c_str()), "Negative daughter E_{clus}/p;E/p (arb.);Counts", 50, 0, 2);
  h_eop_p = new TH1F(Form("%seop_p", histo_prefix.c_str()), "Positive daughter E_{clus}/p;E/p (arb.);Counts", 50, 0, 2);
  // E_clus
  h_E_clus_m = new TH1F(Form("%sE_clus_m", histo_prefix.c_str()), "Negative daughter E_{cluster};E (GeV);Counts", 50, 0, 10);
  h_E_clus_p = new TH1F(Form("%sE_clus_p", histo_prefix.c_str()), "Positive daughter E_{cluster};E (GeV);Counts", 50, 0, 10);
  // Residuals -- Track-cluster
  h_dEta_trackcluster_m = new TH1F(Form("%sdEta_trackcluster_m", histo_prefix.c_str()), "- daughter d#eta = #eta_{track state}-#eta_{clus};d#eta;Counts", 50, -0.5, 0.5);
  h_dEta_trackcluster_p = new TH1F(Form("%sdEta_trackcluster_p", histo_prefix.c_str()), "+ daughter d#eta = #eta_{track state}-#eta_{clus};d#eta;Counts", 50, -0.5, 0.5);
  h_dPhi_trackcluster_m = new TH1F(Form("%sdPhi_trackcluster_m", histo_prefix.c_str()), "- daughter d#varphi = #varphi_{track state}-#varphi_{clus};d#varphi (rad);Counts", 50, -0.5, 0.5);
  h_dPhi_trackcluster_p = new TH1F(Form("%sdPhi_trackcluster_p", histo_prefix.c_str()), "+ daughter d#varphi = #varphi_{track state}-#varphi_{clus};d#varphi (rad);Counts", 50, -0.5, 0.5);
  h_dZ_trackcluster_m =   new TH1F(Form("%sdZ_trackcluster_m",   histo_prefix.c_str()), "- daughter dZ = Z_{track state}-Z_{clus};dZ (cm);Counts", 80, 10, 10);
  h_dZ_trackcluster_p =   new TH1F(Form("%sdZ_trackcluster_p",   histo_prefix.c_str()), "+ daughter dZ = Z_{track state}-Z_{clus};dZ (cm);Counts", 80, 10, 10);
  // nStates in tracking detectors
  h_nINTT_m_vs_p = new TH2F(Form("%snStatesINTT", histo_prefix.c_str()), "Daughter nStates INTT (- vs +);+ daughter hits;- daughter hits", 5,  -0.5, 5.5,  5,  -0.5, 5.5);
  h_nMVTX_m_vs_p = new TH2F(Form("%snStatesMVTX", histo_prefix.c_str()), "Daughter nStates MVTX (- vs +);+ daughter hits;- daughter hits", 5,  -0.5, 5.5,  5,  -0.5, 5.5);
  h_nTPC_m_vs_p  = new TH2F(Form("%snStatesTPC",  histo_prefix.c_str()), "Daughter nStates TPC (- vs +) ;+ daughter hits;- daughter hits", 55, -0.5, 55.5, 55, -0.5, 55.5);  
  // Angle between daughter tracks
  if(m_cut_on_angle_between_tracks){
    h_dEta_tracks = new TH1F(Form("%sdEta_tracks", histo_prefix.c_str()), "Daughter d#eta = |#eta_{-}-#eta_{+}|;d#eta;Counts", 100, -0.5, 0.5);
    h_dPhi_tracks = new TH1F(Form("%sdPhi_tracks", histo_prefix.c_str()), "Daughter d#varphi = |#varphi_{-}-#varphi_{+}|;d#varphi (rad);Counts", 100, -0.5, 0.5);
  }
  // Mother
  h_mother_mass = new TH1F(Form("%smother_mass", histo_prefix.c_str()), "Mother Mass;Invariant M (GeV/c^{2};Counts)", 50, 0, 10);
  // TH2F
  h_mass_vs_eop = new TH2F(Form("%smothermass_vs_eop", histo_prefix.c_str()), "Invariant mass vs. Daughter E/p;Daughter E/p (arb.);M (GeV/c^{2})", 50, 0, 2, 50, 0, 10);
  h_mass_vs_pt  = new TH2F(Form("%smothermass_vs_pt", histo_prefix.c_str()),  "Invariant mass vs. Daughter p_{T};Daughter p_{T} (GeV/c);M (GeV/c^{2})",  50, 0, 10, 50, 0, 10);

  // Register
  if(use_qa_histomanager) // Default is TRUE 
  {
    hm->registerHisto(h_pt_m);
    hm->registerHisto(h_pt_p);
    hm->registerHisto(h_pt_m_vs_p);
    hm->registerHisto(h_eta_vs_phi_m);
    hm->registerHisto(h_eta_vs_phi_p);
    hm->registerHisto(h_eop_m);
    hm->registerHisto(h_eop_p);
    hm->registerHisto(h_E_clus_m);
    hm->registerHisto(h_E_clus_p);
    hm->registerHisto(h_dEta_trackcluster_m);
    hm->registerHisto(h_dEta_trackcluster_p);
    hm->registerHisto(h_dPhi_trackcluster_m);
    hm->registerHisto(h_dPhi_trackcluster_p);
    hm->registerHisto(h_dZ_trackcluster_m);
    hm->registerHisto(h_dZ_trackcluster_p);
    hm->registerHisto(h_nINTT_m_vs_p);
    hm->registerHisto(h_nMVTX_m_vs_p);    
    hm->registerHisto(h_nTPC_m_vs_p);    
    hm->registerHisto(h_mother_mass);
    hm->registerHisto(h_mass_vs_eop);
    hm->registerHisto(h_mass_vs_pt);
    if(m_cut_on_angle_between_tracks)
    {
    hm->registerHisto(h_dEta_tracks);
    hm->registerHisto(h_dPhi_tracks);
    }
  }
}

void QAKFParticleTracktoCalo::FillHistos(
  std::vector<QAKFParticleTracktoCalo::MatchedCluster*> v_matched_clusterstracks,
  float mother_mass,
  std::pair<float,float> anglesbetweendaughters
){
  float pt_p = -1;
  float pt_m = -1;
  unsigned int nStatesINTT_p = 5;
  unsigned int nStatesMVTX_p = 5; 
  unsigned int nStatesTPC_p  = 60; 
  unsigned int nStatesINTT_m = 5; 
  unsigned int nStatesMVTX_m = 5; 
  unsigned int nStatesTPC_m  = 60;
  for(const auto& MaCl : v_matched_clusterstracks){
    //Positive daughters
    if(MaCl->track_charge>0){
      // std::cout << "Filling histos for positrons..." << std::endl;
      pt_p = MaCl->track_pt_svtx;
      h_pt_p->Fill(pt_p);
      float eop = (MaCl->cluster_E/MaCl->track_p_svtx);
      h_eop_p->Fill(eop);
      h_eta_vs_phi_p->Fill(MaCl->track_phi_svtx, MaCl->track_eta_svtx);
      h_E_clus_p->Fill(MaCl->cluster_E);
      h_dEta_trackcluster_p->Fill(MaCl->deta_trackcluster);
      h_dPhi_trackcluster_p->Fill(MaCl->dphi_trackcluster);
      h_dZ_trackcluster_p->Fill(MaCl->dz_trackcluster);
      h_mass_vs_eop->Fill(eop, mother_mass);
      h_mass_vs_pt->Fill(pt_p, mother_mass);
      // nStates
      nStatesINTT_p = MaCl->nStatesINTT;
      nStatesMVTX_p = MaCl->nStatesMVTX;
      nStatesTPC_p =  MaCl->nStatesTPC;
    }
    //Negative daughters
    if(MaCl->track_charge<0){
      // std::cout << "Filling histos for electrons..." << std::endl;
      pt_m = MaCl->track_pt_svtx;
      h_pt_m->Fill(pt_p);
      float eop = (MaCl->cluster_E/MaCl->track_p_svtx);
      h_eop_m->Fill(eop);
      h_eta_vs_phi_m->Fill(MaCl->track_phi_svtx, MaCl->track_eta_svtx);
      h_E_clus_m->Fill(MaCl->cluster_E);
      h_dEta_trackcluster_m->Fill(MaCl->deta_trackcluster);
      h_dPhi_trackcluster_m->Fill(MaCl->dphi_trackcluster);
      h_dZ_trackcluster_m->Fill(MaCl->dz_trackcluster);
      h_mass_vs_eop->Fill(eop, mother_mass);
      h_mass_vs_pt->Fill(pt_m, mother_mass);
      // nStates
      nStatesINTT_m = MaCl->nStatesINTT;
      nStatesMVTX_m = MaCl->nStatesMVTX;
      nStatesTPC_m =  MaCl->nStatesTPC;
    }
  }

  // pT correlation -- both daughters have to have matches in the EMCal 
  if(v_matched_clusterstracks.size()==2){
    h_pt_m_vs_p->Fill(pt_p, pt_m); // x axis +, y axis -
    h_nINTT_m_vs_p->Fill(nStatesINTT_p, nStatesINTT_m);
    h_nMVTX_m_vs_p->Fill(nStatesMVTX_p, nStatesMVTX_m);
    h_nTPC_m_vs_p->Fill(nStatesTPC_p, nStatesTPC_m);
  }
  // Fill histos that only need to be filled ONCE, even if multiple daughters get matched to clusters
  h_mother_mass->Fill(mother_mass); //Need to double check this will work
  if(m_cut_on_angle_between_tracks){
    h_dEta_tracks->Fill(anglesbetweendaughters.first);
    h_dPhi_tracks->Fill(anglesbetweendaughters.second);
  }
}

void QAKFParticleTracktoCalo::WriteHistos(){ // This only runs if use_qa_histomanager is set to false
  m_outfile->cd();

  // Write histos
  std::cout << "Writing histograms..." << std::endl;
  h_pt_p->Write();
  h_pt_m->Write();
  h_eop_p->Write();
  h_eop_m->Write();
  h_E_clus_p->Write();
  h_E_clus_m->Write();
  h_dEta_trackcluster_p->Write();
  h_dPhi_trackcluster_p->Write();
  h_dZ_trackcluster_p->Write();
  h_dEta_trackcluster_m->Write();
  h_dPhi_trackcluster_m->Write();
  h_dZ_trackcluster_m->Write();
  h_eta_vs_phi_p->Write();
  h_eta_vs_phi_m->Write();
  h_pt_m_vs_p->Write();
  h_nINTT_m_vs_p->Write();
  h_nMVTX_m_vs_p->Write();
  h_nTPC_m_vs_p->Write();
  // Mother, correlation
  h_mother_mass->Write();
  h_mass_vs_eop->Write();
  h_mass_vs_pt->Write();
  // Optional
  if(m_cut_on_angle_between_tracks){
    h_dEta_tracks->Write();
    h_dPhi_tracks->Write();
  }
}


//--------------------------------------------------------------------------------
//-----------------------------HELPER FUNCTIONS-----------------------------------
//--------------------------------------------------------------------------------


// ______________________________ MATCHING DONE HERE ________________________________________
QAKFParticleTracktoCalo::MatchedCluster* QAKFParticleTracktoCalo::GetMatchedCluster(SvtxTrackState *thisState, SvtxVertex *m_vertex){
  
  // Objects
  MatchedCluster* MaCl = new MatchedCluster;
  RawCluster *cluster;
  bool is_match;
  int index;

  // Make sure there's a state
  if (thisState == nullptr || m_vertex==nullptr) return MaCl;

  // Vertex
  float vx = m_vertex->get_x();
  float vy = m_vertex->get_y();
  float vz = m_vertex->get_z(); 
  
  // Tracking variables
  float _track_phi_emc = std::numeric_limits<float>::quiet_NaN();
  float _track_eta_emc = std::numeric_limits<float>::quiet_NaN();
  float _track_x_emc = std::numeric_limits<float>::quiet_NaN();
  float _track_y_emc = std::numeric_limits<float>::quiet_NaN();
  float _track_z_emc = std::numeric_limits<float>::quiet_NaN();
  
  // EMCal vectors
  std::vector<float> v_emcal_phi;
  std::vector<float> v_emcal_eta;
  std::vector<float> v_emcal_clusE;
  std::vector<float> v_emcal_dphi;
  std::vector<float> v_emcal_deta;
  std::vector<float> v_emcal_dr;
  std::vector<float> v_emcal_dz;
  
  //Track state kinematics
  _track_x_emc = thisState->get_x()-vx;
  _track_y_emc = thisState->get_y()-vy;
  _track_z_emc = thisState->get_z()-vz;
  _track_phi_emc = std::atan2(_track_y_emc, _track_x_emc);
  _track_eta_emc = std::asinh(_track_z_emc / std::sqrt((_track_x_emc * _track_x_emc) + (_track_y_emc * _track_y_emc)));
 
  // Set variables for matching
  is_match = false;
  index = -1;
  // Create objects, containers, iterators for clusters
  cluster = nullptr;
  RawClusterContainer::Range begin_end_EMC = clustersEM->getClusters();
  RawClusterContainer::Iterator clusIter_EMC;

  // Loop over the EMCal clusters
  for (clusIter_EMC = begin_end_EMC.first; clusIter_EMC != begin_end_EMC.second; ++clusIter_EMC)
  {
    //Get cluster with respect to PV, not (0,0,0)
    cluster = clusIter_EMC->second;
    CLHEP::Hep3Vector v_vertex;
    v_vertex.set(vx,vy,vz);
    CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetEVec(*cluster, v_vertex);

    bool skipthis = false;
    
    //CUT -- CLUSTER ENERGY
    float cluster_energy = E_vec_cluster.mag();
    if (cluster_energy < m_emcal_e_low_cut){
      skipthis=true; failedcut_energy++; 
      std::cout << "Bad energy: " << cluster_energy << std::endl;
    } //continue;}

    //CUT -- DPHI 
    float clPhi = E_vec_cluster.phi();
    float dphi = PiRange(_track_phi_emc - clPhi); 
    if (dphi > m_dphi_cut_high || dphi < m_dphi_cut_low){
      skipthis=true; failedcut_dPhi++; 
      std::cout << "Bad dPhi:   " << dphi << std::endl;
    } // continue;}

    //CUT -- DZ
    float clZ = cluster->get_z() - vz;
    float dz = _track_z_emc - clZ;
    if (dz > m_dz_cut_high || dz < m_dz_cut_low){
      skipthis=true; failedcut_dZ++; 
      std::cout << "Bad dZ:     " << dz << std::endl;
     }// continue;}
   
    // CUT -- DETA    // Currently not being used, but want to include functionality now 
    float clEta = E_vec_cluster.pseudoRapidity();
    float deta = _track_eta_emc - clEta; 
    // if (deta > m_deta_cut_high || deta < m_deta_cut_low) continue;
    
    if(skipthis) continue;

    //Calculate dr
    float tmparg = caloRadiusEMCal * dphi;
    float dr = std::sqrt((tmparg * tmparg) + (dz * dz));  // sqrt((R*dphi)^2 + (dz)^2
    // float dr = sqrt((dphi*dphi + deta*deta)); //previous version
    
    // Add potential match's information to vectors
    v_emcal_phi.push_back(clPhi);
    v_emcal_eta.push_back(clEta);
    v_emcal_clusE.push_back(cluster_energy);
    v_emcal_dphi.push_back(dphi);
    v_emcal_dz.push_back(dz);
    v_emcal_deta.push_back(deta);
    v_emcal_dr.push_back(dr);
    is_match = true;
    number_of_matches++;
  }

  // Find the closest match from all potential matches
  if (is_match == true)
  {
    float tmp = 99999;
    for (long unsigned int i = 0; i < v_emcal_dr.size(); i++)
    {
      if (v_emcal_dr[i] < tmp)
      {
        index = i;
        tmp = v_emcal_dr[i];
      }
    }
  }

  if (is_match == false)
  {
    return MaCl;
  }
  else
  {
    MaCl->track_phi_emc=_track_phi_emc;
    MaCl->track_eta_emc=_track_eta_emc;
    MaCl->track_x_emc=_track_x_emc;
    MaCl->track_y_emc=_track_y_emc;
    MaCl->track_z_emc=_track_z_emc;
    MaCl->cluster_phi=v_emcal_phi[index];
    MaCl->cluster_eta=v_emcal_eta[index];
    MaCl->cluster_E=v_emcal_clusE[index];
    MaCl->dphi_trackcluster=v_emcal_dphi[index];
    MaCl->deta_trackcluster=v_emcal_deta[index];
    MaCl->dz_trackcluster=v_emcal_dz[index];
    MaCl->matched=true;
    return MaCl;
  }
}

// ______________________________ GRABS ALL DAUGHTERS IN EVENT i ________________________________________
std::vector<std::pair<SvtxTrack*,KFParticle*>> QAKFParticleTracktoCalo::GetVectorOfDaughters(KFParticle_Container *thiskfpContainer, SvtxTrackMap *thisTrackMap, std::vector<int> track_ids){

  std::vector<std::pair<SvtxTrack*,KFParticle*>> v_out = {};

  for(const auto& daughterID : track_ids){
    for (auto &iter : *thiskfpContainer){
      KFParticle *thisdaughter = iter.second;
      int thisID = thisdaughter->Id();

      if(thisID == daughterID){ 
        // Get the corresponding track
        SvtxTrack *thisTrack = nullptr;
        auto it_track = thisTrackMap->find(thisID);
        if (it_track != thisTrackMap->end()) thisTrack = it_track->second;
        
        // Add Track and KFP object to the vector
        if(thisTrack!=nullptr && thisdaughter!=nullptr) v_out.push_back({thisTrack,thisdaughter});
        break;
      }
    }
  }

  return v_out;
}

// ______________________________ FOR CUT ON ANGLE BETWEEN DAUGHTERS ________________________________________
std::pair<float,float> QAKFParticleTracktoCalo::GetEtaAndPhiBetweenDaughters(std::vector<std::pair<SvtxTrack*,KFParticle*>> v_daughters){
  
  std::pair<float,float> p_out = {std::numeric_limits<float>::quiet_NaN(),std::numeric_limits<float>::quiet_NaN()};
  if(v_daughters.size() != 2){return p_out;}

  KFParticle *positron; 
  KFParticle *electron;
  float dEta_tracks=9999;
  float dPhi_tracks=9999;
  
  // Make sure they are opposite charges
  if(v_daughters[0].first->get_charge() + v_daughters[1].first->get_charge() != 0){
    return p_out;
  }
  
  // Store the daughters as positrons or electrons based on their charge
  if(v_daughters[0].first->get_charge()>0)
  {
    positron = v_daughters[0].second;
    electron = v_daughters[1].second;
  }
  else
  {
    positron = v_daughters[1].second;
    electron = v_daughters[0].second;
  }
  
  //Get the phi and eta values
  float kfpstphi_p, kfpstphi_m, kfpsteta_p, kfpsteta_m, TempError;
  positron->GetPhi(kfpstphi_p, TempError);
  electron->GetPhi(kfpstphi_m, TempError);
  positron->GetEta(kfpsteta_p, TempError);
  electron->GetEta(kfpsteta_m, TempError);

  // Get dEta and dPhi
  dEta_tracks = kfpsteta_m-kfpsteta_p;
  dPhi_tracks = PiRange(kfpstphi_m-kfpstphi_p);
  p_out = {dEta_tracks,dPhi_tracks};

  // Return statement
  return p_out;
}

bool QAKFParticleTracktoCalo::GoodAngleBetweenDaughters(std::pair<float,float> dEtadPhi){

  // Check Eta and Phi
  float dEta_tracks = dEtadPhi.first;
  float dPhi_tracks = dEtadPhi.second;
  
  if((float) std::abs(dEta_tracks) > m_dEta_Tracks) {return false;}
  if((float) std::abs(dPhi_tracks) > m_dPhi_Tracks) {return false;}

  // Return true if the daughters pass both checks 
  return true;
}


// ______________________________ FOR CUT ON MAX MVTX HITS ___________________________________________
// bool QAKFParticleTracktoCalo::HitsMVTXisGood(PHCompositeNode *topNode, const SvtxTrack *track){
//   int nStatesMVTX = 0;
//   int nStatesINTT = 0;
//   int nStatesTPC  = 0;
//   int nStatesMVTX = GetDetectorHits(topNode, track, nStatesMVTX, nStatesINTT, nStatesTPC);
//   if(m_verbosity > 5){
//     std::cout << "Number of MVTX hits: " << nStatesMVTX << std::endl; 
//     std::cout << "Number of INTT hits: " << nStatesINTT << std::endl; 
//     std::cout << "Number of TPC hits : " << nStatesTPC  << std::endl; 

//   }
//   if(nStatesMVTX > m_MVTX_maxhits )  return false;
//   else if(nStatesMVTX <= m_MVTX_maxhits ) return true;
//   else return false;
// }

void QAKFParticleTracktoCalo::GetDetectorHits(
  const SvtxTrack *track, 
  unsigned int &nStatesMVTX, 
  unsigned int &nStatesINTT,
  unsigned int &nStatesTPC
){
  // dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name_nTuple);
  // if (!dst_trackmap)
  // {
  //   std::cout << "KFParticle truth matching: " << m_trk_map_node_name_nTuple << " does not exist" << std::endl;
  // }

  // detector_nStates_MVTX[daughter_id] = 0;
  // detector_nStates_INTT[daughter_id] = 0;
  // detector_nStates_TPC[daughter_id] = 0;
  // detector_nStates_TPOT[daughter_id] = 0;

  // TrackSeed *silseed = track->get_silicon_seed();
  // // TrackSeed *tpcseed = track->get_tpc_seed();

  // if (!silseed){ return -1;}
  // else
  // {
  //   for (auto cluster_iter = silseed->begin_cluster_keys(); cluster_iter != silseed->end_cluster_keys(); ++cluster_iter)
  //   {
  //     const auto &cluster_key = *cluster_iter;
  //     const auto trackerID = TrkrDefs::getTrkrId(cluster_key);

  //     detector_layer[daughter_id].push_back(TrkrDefs::getLayer(cluster_key));

  //     unsigned int staveId;
  //     unsigned int chipId;
  //     unsigned int ladderZId;
  //     unsigned int ladderPhiId;
  //     unsigned int sectorId;
  //     unsigned int side;
  //     staveId = chipId = ladderZId = ladderPhiId = sectorId = side = std::numeric_limits<unsigned int>::quiet_NaN();

  //     if (trackerID == TrkrDefs::mvtxId)
  //     {
  //       staveId = MvtxDefs::getStaveId(cluster_key);
  //       chipId = MvtxDefs::getChipId(cluster_key);
  //       ++detector_nStates_MVTX[daughter_id];
  //     }
  //     else if (Use["INTT"] && trackerID == TrkrDefs::inttId)
  //     {
  //       ladderZId = InttDefs::getLadderZId(cluster_key);
  //       ladderPhiId = InttDefs::getLadderPhiId(cluster_key);
  //       ++detector_nStates_INTT[daughter_id];
  //     }

  //     if (m_get_detailed_tracking)
  //     {
  //       mvtx_staveID[daughter_id].push_back(staveId);
  //       mvtx_chipID[daughter_id].push_back(chipId);
  //       intt_ladderZID[daughter_id].push_back(ladderZId);
  //       intt_ladderPhiID[daughter_id].push_back(ladderPhiId);
  //       tpc_sectorID[daughter_id].push_back(sectorId);
  //       tpc_side[daughter_id].push_back(side);
  //     }
  //   }
  // }

  // return detector_nStates_MVTX;

  // if (tpcseed)
  // {
  //   for (auto cluster_iter = tpcseed->begin_cluster_keys(); cluster_iter != tpcseed->end_cluster_keys(); ++cluster_iter)
  //   {
  //     const auto &cluster_key = *cluster_iter;
  //     const auto trackerID = TrkrDefs::getTrkrId(cluster_key);

  //     detector_layer[daughter_id].push_back(TrkrDefs::getLayer(cluster_key));

  //     unsigned int staveId;
  //     unsigned int chipId;
  //     unsigned int ladderZId;
  //     unsigned int ladderPhiId;
  //     unsigned int sectorId;
  //     unsigned int side;
  //     staveId = chipId = ladderZId = ladderPhiId = sectorId = side = std::numeric_limits<unsigned int>::quiet_NaN();

  //     if (Use["TPC"] && trackerID == TrkrDefs::tpcId)
  //     {
  //       sectorId = TpcDefs::getSectorId(cluster_key);
  //       side = TpcDefs::getSide(cluster_key);
  //       ++detector_nStates_TPC[daughter_id];
  //     }
  //     else if (Use["TPOT"] && trackerID == TrkrDefs::micromegasId)
  //     {
  //       ++detector_nStates_TPOT[daughter_id];
  //     }

  //     if (m_get_detailed_tracking)
  //     {
  //       mvtx_staveID[daughter_id].push_back(staveId);
  //       mvtx_chipID[daughter_id].push_back(chipId);
  //       intt_ladderZID[daughter_id].push_back(ladderZId);
  //       intt_ladderPhiID[daughter_id].push_back(ladderPhiId);
  //       tpc_sectorID[daughter_id].push_back(sectorId);
  //       tpc_side[daughter_id].push_back(side);
  //     }
  //   }
  // }

  for (auto state_iter = track->begin_states();
       state_iter != track->end_states();
       ++state_iter
  ){

    SvtxTrackState *tstate = state_iter->second;
    if (tstate->get_pathlength() != 0)  // The first track state is an extrapolation so has no cluster
    {
      auto stateckey = tstate->get_cluskey();
      TrkrCluster *cluster = m_dst_clustermap->findCluster(stateckey);
      if (!cluster) continue; // do not have associated cluster, could be track states projected to calo system 
      // auto global = geometry->getGlobalPosition(stateckey, cluster);

      // if (m_get_detailed_tracking)
      // {
      //   residual_x[daughter_id].push_back(global.x() - tstate->get_x());
      //   residual_y[daughter_id].push_back(global.y() - tstate->get_y());
      //   residual_z[daughter_id].push_back(global.z() - tstate->get_z());
      // }

      uint8_t id = TrkrDefs::getTrkrId(stateckey);

      switch (id)
      {
      case TrkrDefs::mvtxId:
        // ++detector_nStates_MVTX[daughter_id];
        nStatesMVTX++;
        break;
      case TrkrDefs::inttId:
        // ++detector_nStates_INTT[daughter_id];
        nStatesINTT++;
        break;
      case TrkrDefs::tpcId:
        // ++detector_nStates_TPC[daughter_id];
        nStatesTPC++;
        break;
      // case TrkrDefs::micromegasId:
      //   // ++detector_nStates_TPOT[daughter_id];
      //   break;
      default:
        // std::cout << "Cluster key doesnt match a tracking system, could be related with projected track state to calorimeter system" << std::endl;
        break;
      }
    }
  }

}

//____________________________________ MISC. FUNCTIONS ______________________________________
std::string QAKFParticleTracktoCalo::GetOutputFileNameAndPath(std::string string_in)
{
  std::string no_ext = string_in;
  const std::string ext = ".root";
  if (string_in.size() >= ext.size() && string_in.compare(string_in.size() - ext.size(), ext.size(), ext) == 0){
    no_ext = string_in.substr(0, string_in.size() - ext.size());
  }

  std::string filename = no_ext;
  std::string path = "";
  std::size_t filestart = no_ext.find_last_of("/");
  if(filestart!=std::string::npos){
    filename = no_ext.substr(filestart+1);
    path = no_ext.substr(0,filestart+1);
  }

  m_outputpath = path;
  m_outtrailer = "_" + m_runnumber + ext;
  std::string file = m_outputpath + filename + m_outtrailer;
  
  // std::cout << "path:    " << m_outputpath << std::endl;
  // std::cout << "name:    " << filename << std::endl;
  // std::cout << "trailer: " << m_outtrailer << std::endl;
  std::cout << "file out: " << file << std::endl;
  
  return file;
}
