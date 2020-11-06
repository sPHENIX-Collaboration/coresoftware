#ifndef KFParticle_truthAndDetTools_H__
#define KFParticle_truthAndDetTools_H__

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4eval/SvtxTrackEval.h>
#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxTruthEval.h>
#include <g4eval/SvtxVertexEval.h>
#include <g4eval/SvtxEvalStack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxVertex.h>

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

    KFParticle_truthAndDetTools(); //Constructor

    ~KFParticle_truthAndDetTools(); //Destructor

    SvtxTrack* getTrack( unsigned int track_id, SvtxTrackMap *trackmap );
    void initializeTruthBranches( TTree *m_tree,int daughter_id );
    void fillTruthBranch( PHCompositeNode *topNode, TTree *m_tree, KFParticle daughter, int daughter_id );

    void initializeDetectorBranches( TTree *m_tree, int daughter_id );
    void initializeSubDetectorBranches( TTree *m_tree, std::string detectorName, int daughter_id );
    void fillDetectorBranch( PHCompositeNode *topNode, TTree *m_tree, KFParticle daughter, int daughter_id );

  private:

    TTree *m_tree;

    float m_true_daughter_vertex_x[20];
    float m_true_daughter_vertex_y[20];
    float m_true_daughter_vertex_z[20];
    float m_true_daughter_px[20];
    float m_true_daughter_py[20];
    float m_true_daughter_pz[20];
    float m_true_daughter_p[20];
    float m_true_daughter_pt[20];
    int m_true_daughter_id[20];

    std::vector<float> detector_local_x[20]; // 7 subdetector including outer and inner hcal plus 4th tracker
    std::vector<float> detector_local_y[20];
    std::vector<float> detector_local_z[20];
    std::vector<int> detector_layer[20];
    std::vector<int> mvtx_staveID[20];
    std::vector<int> mvtx_chipID[20];
    std::vector<int> intt_ladderZID[20];
    std::vector<int> intt_ladderPhiID[20];
    std::vector<int> tpc_sectorID[20];
    std::vector<int> tpc_side[20];

 protected:

    SvtxEvalStack *m_svtx_evalstack;
    SvtxTrackEval *trackeval;
    SvtxClusterEval *clustereval;
    SvtxTruthEval *trutheval;
    SvtxVertexEval *vertexeval;

    SvtxTrackMap *dst_trackmap;
    SvtxTrack *track;

    PHG4Particle* g4particle;
    PHG4VtxPoint* g4vertex_point;
    SvtxVertex *vertex;

    TrkrClusterContainer* dst_clustermap;
};


#endif
