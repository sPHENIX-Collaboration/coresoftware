#ifndef QAKFPARTICLETRACKTOCALO_H
#define QAKFPARTICLETRACKTOCALO_H

#include <fun4all/SubsysReco.h>
// Calo
#include <calobase/RawTowerDefs.h>
#include <calobase/RawClusterDefs.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawClusterUtility.h> // For clusters
// Tower stuff
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainerv4.h>
#include <calobase/TowerInfov4.h>
#include <calobase/TowerInfoDefs.h>
// Vertex
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>
// KFP
#include <KFParticle.h>
#include <kfparticle_sphenix/KFParticle_Container.h>
#include <kfparticle_sphenix/KFParticle_Tools.h>
// Tracking
#include <trackbase_historic/SvtxTrack.h>  // for SvtxTrack
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/InttDefs.h>              // for getLadderPhiId, getLad...
#include <trackbase/MvtxDefs.h>              // for getChipId, getStaveId
#include <trackbase/TpcDefs.h>               // for getSectorId, getSide
#include <trackbase/TrkrCluster.h>           // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>  // for TrkrClusterContainer
#include <trackbase/TrkrDefs.h>              // for getLayer, getTrkrId
// ROOT things
#include <TH1F.h>
#include <TH2F.h>
// Misc.
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cmath>

class PHCompositeNode;
class Fun4AllHistoManager;
class SvtxTrackMap;
class SvtxTrack;
class SvtxVertexMap;
class SvtxVertex;
class TH1;
class TH2;
class TFile;
class QAKFParticleTracktoCalo : public SubsysReco
{
 public:

 // CONSTRUCTOR
  QAKFParticleTracktoCalo(const std::string &name = "QAKFParticleTracktoCalo", const std::string node_name = "");
  ~QAKFParticleTracktoCalo() override = default;

 // OUTPUT FILE -- histo manager is default, but a custom outfile can be used if use_qa_histomanager=false
  void set_histo_prefix(std::string in){histo_prefix = in;}
  void setOutfileName(std::string name_in){
    if(name_in.substr(name_in.size() - 5) != ".root") m_outfilename = name_in.substr(0, name_in.size() - 5);
    else m_outfilename = name_in;
  }
  void dontUseHistoManager(){use_qa_histomanager=false;}
  
 // SETTING CUTS
  void setTrackClusterDeltaZCut(float in){ m_dz_cut_low = in; m_dz_cut_high = in;}
  void setTrackClusterDeltaZMin(float in){ m_dz_cut_high = in;}
  void setTrackClusterDeltaZMax(float in){ m_dz_cut_low = in;}
  void setTrackClusterDeltaPhiCut(float in){ m_dphi_cut_low = in; m_dphi_cut_high = in;}
  void setTrackClusterDeltaPhiMin(float in){ m_dphi_cut_high = in;}
  void setTrackClusterDeltaPhiMax(float in){ m_dphi_cut_low = in;}
  void setClusterEnergyMin(float in){m_emcal_e_low_cut = in;}
  void setProjectionRadius(float in){caloRadiusEMCal = in;}
  // OPT-IN CUTS -- OFF BY DEFAULT -- ONLY ABLE TO BE APPLIED IF THERE ARE TWO DAUGHTER TRACKS
  void setCutOnEtaBetweenTracks(float in){
    m_cut_on_angle_between_tracks = true;
    m_dEta_Tracks=in;
  }
  void setCutOnPhiBetweenTracks(float in){
    m_cut_on_angle_between_tracks = true;
    m_dPhi_Tracks=in;
  }
  void setMaximumHitsMVTX(int in){ 
    m_cut_on_MaxHitsMVTX = true;
    m_MVTX_maxhits = in;
  }

 // VERBOSITY 
  int m_verbosity{0};
  void setVerbosity(int verb) { m_verbosity = verb; }

private:
 
 // STRUCT FOR MATCHING
  struct MatchedCluster{
    // Mother info
    float mother_mass{std::numeric_limits<float>::quiet_NaN()};
    // Info from track, not projected track state
    float track_charge{std::numeric_limits<float>::quiet_NaN()};
    float track_pt_svtx{std::numeric_limits<float>::quiet_NaN()};
    float track_p_svtx{std::numeric_limits<float>::quiet_NaN()};
    float track_eta_svtx{std::numeric_limits<float>::quiet_NaN()};
    float track_phi_svtx{std::numeric_limits<float>::quiet_NaN()};
    float track_phi_KFP{std::numeric_limits<float>::quiet_NaN()};
    float track_eta_KFP{std::numeric_limits<float>::quiet_NaN()};
    unsigned int nStatesINTT={std::numeric_limits<unsigned int>::quiet_NaN()};
    unsigned int nStatesMVTX={std::numeric_limits<unsigned int>::quiet_NaN()};
    unsigned int nStatesTPC={std::numeric_limits<unsigned int>::quiet_NaN()};
    // Info from projected track state and 
    float track_phi_emc{std::numeric_limits<float>::quiet_NaN()};
    float track_eta_emc{std::numeric_limits<float>::quiet_NaN()};
    float track_x_emc{std::numeric_limits<float>::quiet_NaN()};
    float track_y_emc{std::numeric_limits<float>::quiet_NaN()};
    float track_z_emc{std::numeric_limits<float>::quiet_NaN()};
    float cluster_eta{std::numeric_limits<float>::quiet_NaN()};
    float cluster_phi{std::numeric_limits<float>::quiet_NaN()};
    float cluster_E{std::numeric_limits<float>::quiet_NaN()};
    float dphi_trackcluster{std::numeric_limits<float>::quiet_NaN()};
    float deta_trackcluster{std::numeric_limits<float>::quiet_NaN()};
    float dz_trackcluster{std::numeric_limits<float>::quiet_NaN()};
    bool matched{false};
  };

 // OUTPUT FILE
  bool use_qa_histomanager{true};
  std::string histo_prefix = "h_tracktocaloKFPQA_";
  TFile *m_outfile{nullptr};
  std::string m_outputpath; 
  std::string m_outfilename{"track_calo_matching_QA"};
  std::string m_outtrailer{".root"};

 // HISTOS
  // Negative daughters
  TH1F *h_pt_m{nullptr};
  TH1F *h_eop_m{nullptr};
  TH1F *h_E_clus_m{nullptr};
  TH1F *h_dEta_trackcluster_m{nullptr};
  TH1F *h_dPhi_trackcluster_m{nullptr};
  TH1F *h_dZ_trackcluster_m{nullptr};
  TH2F *h_eta_vs_phi_m{nullptr};
  // Positive daughters
  TH1F *h_pt_p{nullptr};
  TH1F *h_eop_p{nullptr};
  TH1F *h_E_clus_p{nullptr};
  TH1F *h_dEta_trackcluster_p{nullptr};
  TH1F *h_dPhi_trackcluster_p{nullptr};
  TH1F *h_dZ_trackcluster_p{nullptr};
  TH1F *h_nINTT_p{nullptr};
  TH1F *h_nMVTX_p{nullptr};
  TH1F *h_nTPC_p{nullptr};
  // Using info from both tracks
  TH2F *h_eta_vs_phi_p{nullptr};
  TH2F *h_nINTT_m_vs_p{nullptr};
  TH2F *h_nMVTX_m_vs_p{nullptr};
  TH2F *h_nTPC_m_vs_p{nullptr};
  TH1F *h_dEta_tracks{nullptr};
  TH1F *h_dPhi_tracks{nullptr};
  TH2F *h_pt_m_vs_p{nullptr};
  // Mother
  TH1F *h_mother_mass{nullptr};
  // Correlation
  TH2F *h_mass_vs_eop{nullptr};
  TH2F *h_mass_vs_pt{nullptr};
  
 // FUNCTIONS
  //FUN4ALL 
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void LoadNodes(PHCompositeNode *topNode);
  void EndEvent(PHCompositeNode* topNode);
  // CREATING AND FILLING HISTOS
  void SetUpHistos();
  void WriteHistos();
  void FillHistos(std::vector<QAKFParticleTracktoCalo::MatchedCluster*> v1, float mother_mass, std::pair<float,float> angle);
  // MISC.
  MatchedCluster* GetMatchedCluster(SvtxTrackState *thisState, SvtxVertex *m_vertex);
  std::vector<std::pair<SvtxTrack*,KFParticle*>> GetVectorOfDaughters(KFParticle_Container *kfp, SvtxTrackMap *tm,std::vector<int> ids);
  std::string GetOutputFileNameAndPath(std::string string_in);
  bool GoodAngleBetweenDaughters(std::pair<float,float> dEtadPhi);
  std::pair<float,float> GetEtaAndPhiBetweenDaughters(std::vector<std::pair<SvtxTrack*,KFParticle*>> v);
  int GetHitsMVTX(PHCompositeNode *topNode, const SvtxTrack *track);
  void GetDetectorHits(const SvtxTrack *track, unsigned int &nStatesMVTX, unsigned int &nStatesINTT, unsigned int &nStatesTPC);
  // HELPER
  float PiRange(float deltaPhi){
    if (deltaPhi > M_PI) deltaPhi -= 2 * M_PI;
    if (deltaPhi < -M_PI) deltaPhi += 2 * M_PI;
    return deltaPhi;
  }

 // CUTS
  float caloRadiusEMCal{102.9};
  float m_emcal_e_low_cut{0.2};
  float m_dphi_cut_low{-0.15};
  float m_dphi_cut_high{0.15};
  float m_dz_cut_low{-10};
  float m_dz_cut_high{10};
  bool m_cut_on_angle_between_tracks{false}; //Useful for photon conversion to set to true
  bool m_cut_on_MaxHitsMVTX{false};
  float m_dEta_Tracks{0.04};
  float m_dPhi_Tracks{0.02};
  unsigned int m_MVTX_maxhits{0};

 // OTHER
  int m_runnumber = std::numeric_limits<int>::quiet_NaN();
  void runnumber(const int run) { m_runnumber = run; }
  int m_event{0};
  int number_of_matches{0};

 // CONTAINERS
  RawTowerGeomContainer *EMCalGeo{nullptr};
  RawClusterContainer *clustersEM{nullptr};
  TowerInfoContainer *_towersEM{nullptr};
  KFParticle_Container *m_kfpContainer{nullptr};
  TrkrClusterContainer *m_dst_clustermap{nullptr};
  ActsGeometry *geometry{nullptr};
  SvtxTrackMap *m_trackMap{nullptr};
  SvtxVertexMap *m_vertexMap{nullptr};
  std::string m_trackMapName{"SvtxTrackMap"};
  std::string m_KFParticleNodeName{"reconstructedParticles"};

  // Problem children
  // KFP Container
  int filledKFPcontainer{0};
  int emptyKFPcontainer{0};
  // Tracks
  int failedcut_nMVTX{0};
  int failedcut_dAngle{0};
  int failedcut_vertex{0};
  int failedcut_ndaugh{0};
  int failedcut_nomtch{0};
  // Clusters
  int failedcut_energy{0};
  int failedcut_dPhi{0};
  int failedcut_dZ{0};

};
#endif //QAKFParticleTracktoCalo_H
