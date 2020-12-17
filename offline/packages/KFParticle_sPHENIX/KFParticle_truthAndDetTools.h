#ifndef KFParticle_truthAndDetTools_H__
#define KFParticle_truthAndDetTools_H__

#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTrackEval.h>
#include <g4eval/SvtxTruthEval.h>
#include <g4eval/SvtxVertexEval.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>

using namespace std;

class SvtxEvalStack;
class PHCompositeNode;
class SvtxTrackMap;
class SvtxTrack;
class TrkrClusterContainer;

class TFile;
class TTree;
class KFParticle;
class KFPVertex;

class KFParticle_truthAndDetTools
{
 public:
  KFParticle_truthAndDetTools();  //Constructor

  ~KFParticle_truthAndDetTools();  //Destructor

  SvtxTrack *getTrack(unsigned int track_id, SvtxTrackMap *trackmap);
  void initializeTruthBranches(TTree *m_tree, int daughter_id);
  void fillTruthBranch(PHCompositeNode *topNode, TTree *m_tree, KFParticle daughter, int daughter_id);

  void initializeDetectorBranches(TTree *m_tree, int daughter_id);
  void initializeSubDetectorBranches(TTree *m_tree, string detectorName, int daughter_id);
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

  vector<float> detector_local_x[max_tracks];  // 7 subdetector including outer and inner hcal plus 4th tracker
  vector<float> detector_local_y[max_tracks];
  vector<float> detector_local_z[max_tracks];
  vector<int> detector_layer[max_tracks];
  vector<int> mvtx_staveID[max_tracks];
  vector<int> mvtx_chipID[max_tracks];
  vector<int> intt_ladderZID[max_tracks];
  vector<int> intt_ladderPhiID[max_tracks];
  vector<int> tpc_sectorID[max_tracks];
  vector<int> tpc_side[max_tracks];

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

#endif
