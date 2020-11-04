#ifndef KFParticle_nTuple_H__
#define KFParticle_nTuple_H__

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>

class SvtxEvalStack;
class PHCompositeNode;
class SvtxTrackMap;
class SvtxTrack;
class TrkrClusterContainer;

class TFile;
class TTree;
class KFParticle;
class KFPVertex;

class KFParticle_nTuple
{
  public:

    KFParticle_nTuple(); //Constructor

    ~KFParticle_nTuple(); //Destructor

    void initializeVariables();
    void initializeBranches();
    void fillBranch( PHCompositeNode *topNode,
                     KFParticle motherParticle,
                     KFParticle vertex,
                     std::vector<KFParticle> daughters,
                     std::vector<KFParticle> intermediates,
                     int nPVs, int multiplicity );

    SvtxTrack* getTrack( unsigned int track_id, SvtxTrackMap *trackmap );
    void initializeTruthBranches( int daughter_id );
    void fillTruthBranch( PHCompositeNode *topNode, KFParticle daughter, int daughter_id );

    void initializeDetectorBranches( int daughter_id );
    void initializeSubDetectorBranches( std::string detectorName, int daughter_id );
    void fillDetectorBranch( PHCompositeNode *topNode, KFParticle daughter, int daughter_id );

  private:

    TTree *m_tree;

    float m_calculated_mother_mass;
    float m_calculated_mother_mass_err;
    float m_calculated_mother_decaytime;
    float m_calculated_mother_decaytime_err;
    float m_calculated_mother_decaylength;
    float m_calculated_mother_decaylength_err;
    float m_calculated_mother_dira;
    float m_calculated_mother_fdchi2;
    float m_calculated_mother_ip;
    float m_calculated_mother_ipchi2;
    float m_calculated_mother_x;
    float m_calculated_mother_y;
    float m_calculated_mother_z;
    float m_calculated_mother_px;
    float m_calculated_mother_py;
    float m_calculated_mother_pz;
    float m_calculated_mother_pe;
    float m_calculated_mother_p;
    float m_calculated_mother_p_err;
    float m_calculated_mother_pt;
    float m_calculated_mother_pt_err;
    float m_calculated_mother_s;
    int   m_calculated_mother_q;
    float m_calculated_mother_eta;
    float m_calculated_mother_rapidity;
    float m_calculated_mother_theta;
    float m_calculated_mother_phi;
    float m_calculated_mother_v;
    float m_calculated_mother_chi2;
    int   m_calculated_mother_ndof;
    //float *m_calculated_mother_cov;
    float m_calculated_mother_cov[21];

    float m_calculated_intermediate_mass[8];
    float m_calculated_intermediate_mass_err[8];
    float m_calculated_intermediate_decaytime[8];
    float m_calculated_intermediate_decaytime_err[8];
    float m_calculated_intermediate_decaylength[8];
    float m_calculated_intermediate_decaylength_err[8];
    float m_calculated_intermediate_ip[8];
    float m_calculated_intermediate_ipchi2[8];
    float m_calculated_intermediate_x[8];
    float m_calculated_intermediate_y[8];
    float m_calculated_intermediate_z[8];
    float m_calculated_intermediate_px[8];
    float m_calculated_intermediate_py[8];
    float m_calculated_intermediate_pz[8];
    float m_calculated_intermediate_pe[8];
    float m_calculated_intermediate_p[8];
    float m_calculated_intermediate_p_err[8];
    float m_calculated_intermediate_pt[8];
    float m_calculated_intermediate_pt_err[8];
    float m_calculated_intermediate_s[8];
    float m_calculated_intermediate_q[8];
    float m_calculated_intermediate_eta[8];
    float m_calculated_intermediate_rapidity[8];
    float m_calculated_intermediate_theta[8];
    float m_calculated_intermediate_phi[8];
    float m_calculated_intermediate_v[8];
    float m_calculated_intermediate_chi2[8];
    float m_calculated_intermediate_ndof[8];
    //float *m_calculated_intermediate_cov[8];
    float m_calculated_intermediate_cov[8][21];

    float m_calculated_daughter_mass[20];
    float m_calculated_daughter_ip[20];
    float m_calculated_daughter_ipchi2[20];
    float m_calculated_daughter_x[20];
    float m_calculated_daughter_y[20];
    float m_calculated_daughter_z[20];
    float m_calculated_daughter_px[20];
    float m_calculated_daughter_py[20];
    float m_calculated_daughter_pz[20];
    float m_calculated_daughter_pe[20];
    float m_calculated_daughter_p[20];
    float m_calculated_daughter_pt[20];
    float m_calculated_daughter_s[20];
    int   m_calculated_daughter_q[20];
    float m_calculated_daughter_eta[20];
    float m_calculated_daughter_rapidity[20];
    float m_calculated_daughter_theta[20];
    float m_calculated_daughter_phi[20];
    float m_calculated_daughter_chi2[20];
    int   m_calculated_daughter_ndof[20];
    int   m_calculated_daughter_trid[20];
    //float *m_calculated_daughter_cov[20];
    float m_calculated_daughter_cov[20][21];

    float m_true_daughter_px[20];
    float m_true_daughter_py[20];
    float m_true_daughter_pz[20];

    float m_daughter_dca[99];
    
    float m_calculated_vertex_x;
    float m_calculated_vertex_y;
    float m_calculated_vertex_z;
    float m_calculated_vertex_v;  
    float m_calculated_vertex_chi2;
    float m_calculated_vertex_ndof;
    //float *m_calculated_vertex_cov;
    float m_calculated_vertex_cov[6];
 
    int m_nPVs;
    int m_multiplicity;

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

    bool m_has_intermediates_nTuple;
    int m_num_tracks_nTuple;
    int m_num_intermediate_states_nTuple;
    bool m_truth_matching;
    bool m_detector_info;
    std::string m_mother_name;
    bool m_use_intermediate_name;
    std::string m_intermediate_name_ntuple[99];
    SvtxEvalStack *m_svtx_evalstack;
    TrkrClusterContainer* dst_clustermap;
};


#endif
