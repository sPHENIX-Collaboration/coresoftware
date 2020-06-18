#ifndef KFParticle_nTuple_H__
#define KFParticle_nTuple_H__

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
    void initializeBranches( int nTracks);
    void fillBranch( KFParticle motherParticle,
                     KFParticle vertex,
                     int nTracks,
                     KFParticle daughter_1,
                     KFParticle daughter_2,
                     KFParticle daughter_3,
                     KFParticle daughter_4,
                     int nPVs, int multiplicity );

  private:

    TTree *m_tree;

    float m_calculated_mother_mass;
    float m_calculated_mother_mass_err;
    float m_calculated_mother_lifetime;
    float m_calculated_mother_lifetime_err;
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
    float m_calculated_mother_pt;
    float m_calculated_mother_s;
    int   m_calculated_mother_q;
    float m_calculated_mother_eta;
    float m_calculated_mother_rapidity;
    float m_calculated_mother_theta;
    float m_calculated_mother_phi;
    float m_calculated_mother_chi2;
    int   m_calculated_mother_ndof;
    float *m_calculated_mother_cov;

    float m_calculated_daughter_mass[4];
    float m_calculated_daughter_ip[4];
    float m_calculated_daughter_ipchi2[4];
    float m_calculated_daughter_x[4];
    float m_calculated_daughter_y[4];
    float m_calculated_daughter_z[4];
    float m_calculated_daughter_px[4];
    float m_calculated_daughter_py[4];
    float m_calculated_daughter_pz[4];
    float m_calculated_daughter_pe[4];
    float m_calculated_daughter_p[4];
    float m_calculated_daughter_pt[4];
    float m_calculated_daughter_s[4];
    int   m_calculated_daughter_q[4];
    float m_calculated_daughter_eta[4];
    float m_calculated_daughter_rapidity[4];
    float m_calculated_daughter_theta[4];
    float m_calculated_daughter_phi[4];
    float m_calculated_daughter_chi2[4];
    int   m_calculated_daughter_ndof[4];
    float *m_calculated_daughter_cov[4];
    
    float m_calculated_vertex_x;
    float m_calculated_vertex_y;
    float m_calculated_vertex_z;
    float *m_calculated_vertex_cov;
   
    int m_nPVs;
    int m_multiplicity;
};


#endif
