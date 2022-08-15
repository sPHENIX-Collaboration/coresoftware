#ifndef KFPARTICLESPHENIX_KFPARTICLETRUTHANDDETTOOLS_H
#define KFPARTICLESPHENIX_KFPARTICLETRUTHANDDETTOOLS_H

#include <string>
#include <vector>

class PHCompositeNode;
class PHG4Particle;
class PHG4VtxPoint;
class PHG4TruthInfoContainer;
class PHHepMCGenEvent;
class PHHepMCGenEventMap;
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

namespace HepMC
{
  class GenParticle;
}

class KFParticle_truthAndDetTools
{
 public:
  KFParticle_truthAndDetTools();  //Constructor

  virtual ~KFParticle_truthAndDetTools();  //Destructor

  SvtxTrack *getTrack(unsigned int track_id, SvtxTrackMap *trackmap);
  SvtxVertex *getVertex(unsigned int vertex_id, SvtxVertexMap *vertexmap);
  PHG4Particle *getTruthTrack(SvtxTrack *thisTrack, PHCompositeNode *topNode);

  void initializeTruthBranches(TTree *m_tree, int daughter_id, std::string daughter_number, bool m_constrain_to_vertex_truthMatch);
  void fillTruthBranch(PHCompositeNode *topNode, TTree *m_tree, KFParticle daughter, int daughter_id, KFParticle vertex, bool m_constrain_to_vertex_truthMatch);

  void fillGeant4Branch(PHG4Particle *particle, int daughter_id);
  void fillHepMCBranch(HepMC::GenParticle *particle, int daughter_id);
  int getHepMCInfo(PHCompositeNode *topNode, TTree *m_tree, KFParticle daughter, int daughter_id);

  void initializeCaloBranches(TTree *m_tree, int daughter_id, std::string daughter_number);
  void fillCaloBranch(PHCompositeNode *topNode, TTree *m_tree, KFParticle daughter, int daughter_id);

  void initializeDetectorBranches(TTree *m_tree, int daughter_id, std::string daughter_number);
  void initializeSubDetectorBranches(TTree *m_tree, std::string detectorName, int daughter_id, std::string daughter_number);
  void fillDetectorBranch(PHCompositeNode *topNode, TTree *m_tree, KFParticle daughter, int daughter_id);

  void allPVInfo(PHCompositeNode *topNode, TTree *m_tree,
                 KFParticle motherParticle,
                 std::vector<KFParticle> daughters,
                 std::vector<KFParticle> intermediates);

  void clearVectors();

 protected:
  std::string m_trk_map_node_name_nTuple = "SvtxTrackMap";
  std::string m_vtx_map_node_name_nTuple = "SvtxVertexMap";

  SvtxEvalStack *m_svtx_evalstack = nullptr;
  SvtxClusterEval *clustereval = nullptr;
  SvtxHitEval *hiteval = nullptr;
  SvtxTrackEval *trackeval = nullptr;
  SvtxTruthEval *trutheval = nullptr;
  SvtxVertexEval *vertexeval = nullptr;

  SvtxTrackMap *dst_trackmap = nullptr;
  SvtxTrack *track = nullptr;

  PHG4Particle *g4particle = nullptr;
  PHG4VtxPoint *g4vertex_point = nullptr;

  SvtxVertexMap *dst_vertexmap = nullptr;
  SvtxVertex *vertex = nullptr;

  TrkrClusterContainer *dst_clustermap = nullptr;

  int m_num_tracks_nTuple = 0;
  int m_num_intermediate_states_nTuple = 0;

  static const int max_tracks = 20;

  float m_true_daughter_vertex_x[max_tracks] = {0};
  float m_true_daughter_vertex_y[max_tracks] = {0};
  float m_true_daughter_vertex_z[max_tracks] = {0};
  float m_true_daughter_ip[max_tracks] = {0};
  float m_true_daughter_ip_xy[max_tracks] = {0};
  float m_true_daughter_px[max_tracks] = {0};
  float m_true_daughter_py[max_tracks] = {0};
  float m_true_daughter_pz[max_tracks] = {0};
  float m_true_daughter_p[max_tracks] = {0};
  float m_true_daughter_pt[max_tracks] = {0};
  int m_true_daughter_id[max_tracks] = {0};
  float m_true_daughter_pv_x[max_tracks] = {0};
  float m_true_daughter_pv_y[max_tracks] = {0};
  float m_true_daughter_pv_z[max_tracks] = {0};

  std::vector<int> m_true_daughter_track_history_PDG_ID[max_tracks];
  std::vector<float> m_true_daughter_track_history_PDG_mass[max_tracks];
  std::vector<float> m_true_daughter_track_history_px[max_tracks];
  std::vector<float> m_true_daughter_track_history_py[max_tracks];
  std::vector<float> m_true_daughter_track_history_pz[max_tracks];
  std::vector<float> m_true_daughter_track_history_pE[max_tracks];
  std::vector<float> m_true_daughter_track_history_pT[max_tracks];

  float detector_emcal_deltaphi[max_tracks] = {0};
  float detector_emcal_deltaeta[max_tracks] = {0};
  float detector_emcal_energy_3x3[max_tracks] = {0};
  float detector_emcal_energy_5x5[max_tracks] = {0};
  float detector_emcal_cluster_energy[max_tracks] = {0};
  float detector_ihcal_deltaphi[max_tracks] = {0};
  float detector_ihcal_deltaeta[max_tracks] = {0};
  float detector_ihcal_energy_3x3[max_tracks] = {0};
  float detector_ihcal_energy_5x5[max_tracks] = {0};
  float detector_ihcal_cluster_energy[max_tracks] = {0};
  float detector_ohcal_deltaphi[max_tracks] = {0};
  float detector_ohcal_deltaeta[max_tracks] = {0};
  float detector_ohcal_energy_3x3[max_tracks] = {0};
  float detector_ohcal_energy_5x5[max_tracks] = {0};
  float detector_ohcal_cluster_energy[max_tracks] = {0};

  std::vector<float> detector_local_x[max_tracks];  // 7 subdetector including outer and inner hcal plus 4th tracker
  std::vector<float> detector_local_y[max_tracks];
  std::vector<float> detector_local_z[max_tracks];
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

  PHG4TruthInfoContainer *m_truthinfo = nullptr;
  PHHepMCGenEventMap *m_geneventmap = nullptr;
  PHHepMCGenEvent *m_genevt = nullptr;
};

#endif  //KFPARTICLESPHENIX_KFPARTICLETRUTHANDDETTOOLS_H
