#ifndef KFPARTICLESPHENIX_KFPARTICLETRUTHANDDETTOOLS_H
#define KFPARTICLESPHENIX_KFPARTICLETRUTHANDDETTOOLS_H

#include <KFParticle.h>

#include <string>
#include <vector>

class PHCompositeNode;
class PHG4Particle;
class PHG4VtxPoint;
class SvtxClusterEval;
class SvtxEvalStack;
class SvtxTrackMap;
class SvtxTrack;
class SvtxTruthEval;
class SvtxVertex;
class TrkrClusterContainer;

class TTree;
class KFParticle;

class KFParticle_truthAndDetTools
{
 public:
  KFParticle_truthAndDetTools();  //Constructor

  virtual ~KFParticle_truthAndDetTools();  //Destructor

  SvtxTrack *getTrack(unsigned int track_id, SvtxTrackMap *trackmap);
  void initializeTruthBranches(TTree *m_tree, int daughter_id);
  void fillTruthBranch(PHCompositeNode *topNode, TTree *m_tree, KFParticle daughter, int daughter_id);

  void initializeDetectorBranches(TTree *m_tree, int daughter_id);
  void initializeSubDetectorBranches(TTree *m_tree, std::string detectorName, int daughter_id);
  void fillDetectorBranch(PHCompositeNode *topNode, TTree *m_tree, KFParticle daughter, int daughter_id);

 private:
  TTree *m_tree;

  static const int max_tracks = 20;

  float m_true_daughter_vertex_x[max_tracks] = {0};
  float m_true_daughter_vertex_y[max_tracks] = {0};
  float m_true_daughter_vertex_z[max_tracks] = {0};
  float m_true_daughter_px[max_tracks] = {0};
  float m_true_daughter_py[max_tracks] = {0};
  float m_true_daughter_pz[max_tracks] = {0};
  float m_true_daughter_p[max_tracks] = {0};
  float m_true_daughter_pt[max_tracks] = {0};
  int m_true_daughter_id[max_tracks] = {0};

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

 protected:
  SvtxEvalStack *m_svtx_evalstack = nullptr;
  SvtxClusterEval *clustereval = nullptr;
  SvtxTruthEval *trutheval = nullptr;

  SvtxTrackMap *dst_trackmap = nullptr;
  SvtxTrack *track = nullptr;

  PHG4Particle *g4particle = nullptr;
  PHG4VtxPoint *g4vertex_point = nullptr;
  SvtxVertex *vertex = nullptr;

  TrkrClusterContainer *dst_clustermap = nullptr;
};

#endif //KFPARTICLESPHENIX_KFPARTICLETRUTHANDDETTOOLS_H
