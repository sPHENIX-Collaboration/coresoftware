#ifndef KFPARTICLESPHENIX_KFPARTICLETRUTHANDDETTOOLS_H
#define KFPARTICLESPHENIX_KFPARTICLETRUTHANDDETTOOLS_H

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <trackbase/ActsGeometry.h>

#include <limits>
#include <string>
#include <vector>

class ActsGeometry;
class PHCompositeNode;
class PHG4Particle;
class PHG4VtxPoint;
class PHG4TruthInfoContainer;
class PHHepMCGenEvent;
class PHHepMCGenEventMap;
class MbdVertexMap;
class MbdVertex;
class SvtxClusterEval;
class SvtxEvalStack;
class SvtxHitEval;
class SvtxTrack;
class SvtxTrackEval;
class SvtxTrackMap;
class SvtxTruthEval;
class SvtxVertexMap;
class SvtxVertex;
class SvtxVertexEval;
class TrkrClusterContainer;
class TTree;
class KFParticle;
class GlobalVertex;
class GlobalVertexMap;

namespace HepMC
{
  class GenParticle;
}

class KFParticle_truthAndDetTools
{
 public:
  KFParticle_truthAndDetTools() = default;  // Constructor

  virtual ~KFParticle_truthAndDetTools() = default;

  static SvtxTrack *getTrack(unsigned int track_id, SvtxTrackMap *trackmap);
  static GlobalVertex *getVertex(unsigned int vertex_id, GlobalVertexMap *vertexmap);
  PHG4Particle *getTruthTrack(SvtxTrack *thisTrack, PHCompositeNode *topNode);

  void initializeTruthBranches(TTree *m_tree, int daughter_id, const std::string &daughter_number, bool m_constrain_to_vertex_truthMatch);
  void fillTruthBranch(PHCompositeNode *topNode, TTree *m_tree, const KFParticle &daughter, int daughter_id, const KFParticle &vertex, bool m_constrain_to_vertex_truthMatch);

  void fillGeant4Branch(PHG4Particle *particle, int daughter_id);
  void fillHepMCBranch(HepMC::GenParticle *particle, int daughter_id);
  int getHepMCInfo(PHCompositeNode *topNode, TTree *m_tree, const KFParticle &daughter, int daughter_id);

  void initializeCaloBranches(TTree *m_tree, int daughter_id, const std::string &daughter_number);
  void fillCaloBranch(PHCompositeNode *topNode, TTree *m_tree, const KFParticle &daughter, int daughter_id, bool &isTrackEMCalmatch);

  void initializeDetectorBranches(TTree *m_tree, int daughter_id, const std::string &daughter_number);
  void initializeSubDetectorBranches(TTree *m_tree, const std::string &detectorName, int daughter_id, const std::string &daughter_number);
  void fillDetectorBranch(PHCompositeNode *topNode, TTree *m_tree, const KFParticle &daughter, int daughter_id);

  int getPVID(PHCompositeNode *topNode, const KFParticle &vertex);
  void allPVInfo(PHCompositeNode *topNode, TTree *m_tree,
                 const KFParticle &motherParticle,
                 std::vector<KFParticle> daughters,
                 std::vector<KFParticle> intermediates);

  void clearVectors();

  float get_e3x3(RawCluster *cluster, RawTowerContainer *Towers, int layer); //Nonfunctional; for now, tree fills with NAN
  float get_e5x5(RawCluster *cluster, RawTowerContainer *Towers, int layer);

  float PiRange(float deltaPhi)
  {
    if (deltaPhi > M_PI) deltaPhi -= 2 * M_PI;
    if (deltaPhi < -M_PI) deltaPhi += 2 * M_PI;
    return deltaPhi;
  }

  // Functions to set cuts
  void set_emcal_radius_user(float set_variable) { m_emcal_radius_user = set_variable; }
  void set_emcal_e_low_cut(float set_variable) { m_emcal_e_low_cut = set_variable; }
  void set_dphi_cut_low(float set_variable) { m_dphi_cut_low = set_variable; }
  void set_dphi_cut_high(float set_variable) { m_dphi_cut_high = set_variable; }
  void set_dz_cut_low(float set_variable) { m_dz_cut_low = set_variable; }
  void set_dz_cut_high(float set_variable) { m_dz_cut_high = set_variable; }

 protected:
  bool m_use_mbd_vertex_truth{false};
  bool m_dont_use_global_vertex_truth{false};

  std::string m_trk_map_node_name_nTuple{"SvtxTrackMap"};
  std::string m_vtx_map_node_name_nTuple{"SvtxVertexMap"};

  SvtxEvalStack *m_svtx_evalstack{nullptr};
  SvtxClusterEval *clustereval{nullptr};
  SvtxHitEval *hiteval{nullptr};
  SvtxTrackEval *trackeval{nullptr};
  SvtxTruthEval *trutheval{nullptr};
  SvtxVertexEval *vertexeval{nullptr};

  ActsGeometry *geometry{nullptr};

  SvtxTrackMap *dst_trackmap{nullptr};
  SvtxTrack *track{nullptr};

  PHG4Particle *g4particle{nullptr};
  PHG4VtxPoint *g4vertex_point{nullptr};

  MbdVertexMap *dst_mbdvertexmap{nullptr};

  SvtxVertexMap *dst_vertexmap{nullptr};
  SvtxVertex *vertex{nullptr};

  TrkrClusterContainer *dst_clustermap{nullptr};

  int m_num_tracks_nTuple{0};
  int m_num_intermediate_states_nTuple{0};

  static const int max_tracks{20};

  float m_emcal_radius_user{100.70};
  float m_ihcal_radius_user{117};
  float m_ohcal_radius_user{177.423};

  float m_emcal_e_low_cut{0.2};
  float m_ihcal_e_low_cut{0.01};
  float m_ohcal_e_low_cut{0.01};
  float m_dphi_cut_low{-0.15};
  float m_dphi_cut_high{0.15};
  float m_dz_cut_low{-10};
  float m_dz_cut_high{10};

  float m_true_daughter_vertex_x[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_vertex_y[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_vertex_z[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_ip[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_ip_xy[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_px[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_py[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_pz[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_p[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_pt[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  int m_true_daughter_id[max_tracks]{std::numeric_limits<int>::quiet_NaN()};
  float m_true_daughter_pv_x[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_pv_y[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_daughter_pv_z[max_tracks]{std::numeric_limits<float>::quiet_NaN()};

  std::vector<int> m_true_daughter_track_history_PDG_ID[max_tracks];
  std::vector<float> m_true_daughter_track_history_PDG_mass[max_tracks];
  std::vector<float> m_true_daughter_track_history_px[max_tracks];
  std::vector<float> m_true_daughter_track_history_py[max_tracks];
  std::vector<float> m_true_daughter_track_history_pz[max_tracks];
  std::vector<float> m_true_daughter_track_history_pE[max_tracks];
  std::vector<float> m_true_daughter_track_history_pT[max_tracks];

  float detector_emcal_deltaphi[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_emcal_deltaeta[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_emcal_deltaz[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_emcal_energy_3x3[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_emcal_energy_5x5[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_emcal_cluster_energy[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ihcal_deltaphi[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ihcal_deltaeta[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ihcal_energy_3x3[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ihcal_energy_5x5[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ihcal_cluster_energy[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ohcal_deltaphi[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ohcal_deltaeta[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ohcal_energy_3x3[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ohcal_energy_5x5[max_tracks]{std::numeric_limits<float>::quiet_NaN()};
  float detector_ohcal_cluster_energy[max_tracks]{std::numeric_limits<float>::quiet_NaN()};

  RawTowerGeomContainer *EMCalGeo{nullptr};
  RawClusterContainer *clustersEM{nullptr};
  RawTowerGeomContainer *IHCalGeo{nullptr};
  RawClusterContainer *clustersIH{nullptr};
  RawTowerGeomContainer *OHCalGeo{nullptr};
  RawClusterContainer *clustersOH{nullptr};
  RawTowerContainer *_towersEM{nullptr};
  RawTowerContainer *_towersIH{nullptr};
  RawTowerContainer *_towersOH{nullptr};

  unsigned int detector_nHits_MVTX[max_tracks]{0};
  unsigned int detector_nHits_INTT[max_tracks]{0};
  unsigned int detector_nHits_TPC[max_tracks]{0};
  unsigned int detector_nHits_TPOT[max_tracks]{0};
  unsigned int detector_nStates_MVTX[max_tracks]{0};
  unsigned int detector_nStates_INTT[max_tracks]{0};
  unsigned int detector_nStates_TPC[max_tracks]{0};
  unsigned int detector_nStates_TPOT[max_tracks]{0};
  std::vector<float> residual_x[max_tracks];
  std::vector<float> residual_y[max_tracks];
  std::vector<float> residual_z[max_tracks];
  std::vector<int> detector_layer[max_tracks];
  std::vector<int> mvtx_staveID[max_tracks];
  std::vector<int> mvtx_chipID[max_tracks];
  std::vector<int> intt_ladderZID[max_tracks];
  std::vector<int> intt_ladderPhiID[max_tracks];
  std::vector<int> tpc_sectorID[max_tracks];
  std::vector<int> tpc_side[max_tracks];

  std::vector<float> allPV_x;
  std::vector<float> allPV_y;
  std::vector<float> allPV_z;
  std::vector<float> allPV_mother_IP;
  std::vector<float> allPV_mother_IPchi2;
  std::vector<float> allPV_daughter_IP[max_tracks];
  std::vector<float> allPV_daughter_IPchi2[max_tracks];
  std::vector<float> allPV_intermediates_IP[max_tracks];
  std::vector<float> allPV_intermediates_IPchi2[max_tracks];

  PHG4TruthInfoContainer *m_truthinfo{nullptr};
  PHHepMCGenEventMap *m_geneventmap{nullptr};
  PHHepMCGenEvent *m_genevt{nullptr};
};

#endif  // KFPARTICLESPHENIX_KFPARTICLETRUTHANDDETTOOLS_H
