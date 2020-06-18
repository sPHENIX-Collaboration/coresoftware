#include <KFParticle_nTuple.h>
#include <TTree.h>
#include "KFParticle.h"
#include "KFPVertex.h"

#include <KFParticle_Tools.h>
KFParticle_Tools kfpTupleTools;

KFParticle_nTuple::KFParticle_nTuple(){} //Constructor

KFParticle_nTuple::~KFParticle_nTuple(){} //Destructor

void KFParticle_nTuple::initializeVariables()
{
  m_calculated_mother_mass = -99;
  m_calculated_mother_mass_err = -99;
  m_calculated_mother_lifetime = -99;
  m_calculated_mother_lifetime_err = -99;
  m_calculated_mother_dira = -99;
  m_calculated_mother_fdchi2 = -99;
  m_calculated_mother_ip = -99;
  m_calculated_mother_ipchi2 = -99;
  m_calculated_mother_x = -99;
  m_calculated_mother_y = -99;
  m_calculated_mother_z = -99;
  m_calculated_mother_px = -99;
  m_calculated_mother_py = -99;
  m_calculated_mother_pz = -99;
  m_calculated_mother_pe = -99;
  m_calculated_mother_p = -99;
  m_calculated_mother_pt = -99;
  m_calculated_mother_s = -99;
  m_calculated_mother_q = -99;
  m_calculated_mother_eta = -99;
  m_calculated_mother_rapidity = -99;
  m_calculated_mother_theta = -99;
  m_calculated_mother_phi = -99;
  m_calculated_mother_chi2 = -99;
  m_calculated_mother_ndof = -99;
  //m_calculated_mother_cov = -99;
}

void KFParticle_nTuple::initializeBranches( int nTracks = 2)
{

  m_tree = new TTree("DecayTree", "DecayTree");

  m_tree->Branch( "mother_mass",           &m_calculated_mother_mass,         "mother_mass/F" );
  m_tree->Branch( "mother_massErr",        &m_calculated_mother_mass_err,     "mother_massErr/F" );
  m_tree->Branch( "mother_lifetime",       &m_calculated_mother_lifetime,     "mother_lifetime/F" );
  m_tree->Branch( "mother_lifetimeErr",    &m_calculated_mother_lifetime_err, "mother_lifetimeErr/F" );
  m_tree->Branch( "mother_DIRA",           &m_calculated_mother_dira,         "mother_DIRA/F" );
  m_tree->Branch( "mother_FDchi2",         &m_calculated_mother_fdchi2,       "mother_FDCHI2/F" );
  m_tree->Branch( "mother_IP",             &m_calculated_mother_ip,           "mother_IP/F" );
  m_tree->Branch( "mother_IPCHI2",         &m_calculated_mother_ipchi2,       "mother_IPCHI2/F" );
  m_tree->Branch( "mother_x",              &m_calculated_mother_x,            "mother_x/F" );
  m_tree->Branch( "mother_y",              &m_calculated_mother_y,            "mother_y/F" );
  m_tree->Branch( "mother_z",              &m_calculated_mother_z,            "mother_z/F" );
  m_tree->Branch( "mother_px",             &m_calculated_mother_px,           "mother_px/F" );
  m_tree->Branch( "mother_py",             &m_calculated_mother_py,           "mother_py/F" );
  m_tree->Branch( "mother_pz",             &m_calculated_mother_pz,           "mother_pz/F" );
  m_tree->Branch( "mother_pe",             &m_calculated_mother_pe,           "mother_pe/F" );
  m_tree->Branch( "mother_p",              &m_calculated_mother_p,            "mother_p/F" );
  m_tree->Branch( "mother_pT",             &m_calculated_mother_pt,           "mother_pT/F" );
  m_tree->Branch( "mother_S",              &m_calculated_mother_s,            "mother_S/F" );
  m_tree->Branch( "mother_charge",         &m_calculated_mother_q,            "mother_charge/C" );
  m_tree->Branch( "mother_pseudorapidity", &m_calculated_mother_eta,          "mother_pseudorapidity/F" );
  m_tree->Branch( "mother_rapidity",       &m_calculated_mother_rapidity,     "mother_rapidity/F" );
  m_tree->Branch( "mother_theta",          &m_calculated_mother_theta,        "mother_theta/F" );
  m_tree->Branch( "mother_phi",            &m_calculated_mother_phi,          "mother_phi/F" );
  m_tree->Branch( "mother_chi2",           &m_calculated_mother_chi2,         "mother_chi2/F" );
  m_tree->Branch( "mother_nDoF",           &m_calculated_mother_ndof,         "mother_nDoF/I" );
  m_tree->Branch( "mother_Covariance",     &m_calculated_mother_cov,          "mother_Covariance/F[21]", 21 );

 for (int i = 0; i < nTracks; ++i)
 {
    std::string daughter_number = "daughter_" + std::to_string(i + 1);

    m_tree->Branch( TString(daughter_number) + "_mass",           &m_calculated_daughter_mass[i],     TString(daughter_number) + "_mass/F" );
    m_tree->Branch( TString(daughter_number) + "_IP",             &m_calculated_daughter_ip[i],       TString(daughter_number) + "_IP/F" );
    m_tree->Branch( TString(daughter_number) + "_IPCHI2",         &m_calculated_daughter_ipchi2[i],   TString(daughter_number) + "_IPCHI2/F" );
    m_tree->Branch( TString(daughter_number) + "_x",              &m_calculated_daughter_x[i],        TString(daughter_number) + "_x/F" );
    m_tree->Branch( TString(daughter_number) + "_y",              &m_calculated_daughter_y[i],        TString(daughter_number) + "_y/F" );
    m_tree->Branch( TString(daughter_number) + "_z",              &m_calculated_daughter_z[i],        TString(daughter_number) + "_z/F" );
    m_tree->Branch( TString(daughter_number) + "_px",             &m_calculated_daughter_px[i],       TString(daughter_number) + "_px/F" );
    m_tree->Branch( TString(daughter_number) + "_py",             &m_calculated_daughter_py[i],       TString(daughter_number) + "_py/F" );
    m_tree->Branch( TString(daughter_number) + "_pz",             &m_calculated_daughter_pz[i],       TString(daughter_number) + "_pz/F" );
    m_tree->Branch( TString(daughter_number) + "_pe",             &m_calculated_daughter_pe[i],       TString(daughter_number) + "_pe/F" );
    m_tree->Branch( TString(daughter_number) + "_p",              &m_calculated_daughter_p[i],        TString(daughter_number) + "_p/F" );
    m_tree->Branch( TString(daughter_number) + "_pT",             &m_calculated_daughter_pt[i],       TString(daughter_number) + "_pT/F" );
    m_tree->Branch( TString(daughter_number) + "_S",              &m_calculated_daughter_s[i],        TString(daughter_number) + "_S/F" );
    m_tree->Branch( TString(daughter_number) + "_charge",         &m_calculated_daughter_q[i],        TString(daughter_number) + "_charge/C" );
    m_tree->Branch( TString(daughter_number) + "_pseudorapidity", &m_calculated_daughter_eta[i],      TString(daughter_number) + "_pseudorapidity/F" );
    m_tree->Branch( TString(daughter_number) + "_rapidity",       &m_calculated_daughter_rapidity[i], TString(daughter_number) + "_rapidity/F" );
    m_tree->Branch( TString(daughter_number) + "_theta",          &m_calculated_daughter_theta[i],    TString(daughter_number) + "_theta/F" );
    m_tree->Branch( TString(daughter_number) + "_phi",            &m_calculated_daughter_phi[i],      TString(daughter_number) + "_phi/F" );
    m_tree->Branch( TString(daughter_number) + "_chi2",           &m_calculated_daughter_chi2[i],     TString(daughter_number) + "_chi2/F" );
    m_tree->Branch( TString(daughter_number) + "_nDoF",           &m_calculated_daughter_ndof[i],     TString(daughter_number) + "_nDoF/I" );
    m_tree->Branch( TString(daughter_number) + "_Covariance",     &m_calculated_daughter_cov[i],      TString(daughter_number) + "_Covariance/F[21]", 21 );
}

  m_tree->Branch( "vertex_x",              &m_calculated_vertex_x,            "vertex_x/F" );
  m_tree->Branch( "vertex_y",              &m_calculated_vertex_y,            "vertex_y/F" );
  m_tree->Branch( "vertex_z",              &m_calculated_vertex_z,            "vertex_z/F" );
  m_tree->Branch( "vertex_Covariance",     &m_calculated_vertex_cov,          "vertex_Covariance/F[6]", 6 );

  m_tree->Branch( "nPrimaryVertices",     &m_nPVs,                            "nPrimaryVertices/I" );
  m_tree->Branch( "nEventTracks",         &m_multiplicity,                    "nEventTracks/I" );
}


void KFParticle_nTuple::fillBranch( KFParticle motherParticle, 
                                    KFParticle vertex,
                                    int nTracks, 
                                    KFParticle daughter_1, 
                                    KFParticle daughter_2, 
                                    KFParticle daughter_3, 
                                    KFParticle daughter_4,
                                    int nPVs, int multiplicity )
{

  KFPVertex *kfpVertex = new KFPVertex; 
  kfpVertex->SetXYZ( vertex.Parameters() );
  kfpVertex->SetCovarianceMatrix( vertex.CovarianceMatrix() ); 

  m_calculated_mother_mass         = motherParticle.GetMass();
  m_calculated_mother_mass_err     = motherParticle.GetErrMass();
  m_calculated_mother_lifetime     = motherParticle.GetLifeTime();
  m_calculated_mother_lifetime_err = motherParticle.GetErrLifeTime();
  m_calculated_mother_dira         = kfpTupleTools.eventDIRA( motherParticle, *kfpVertex );
  m_calculated_mother_fdchi2       = kfpTupleTools.flightDistanceChi2( motherParticle, *kfpVertex );
  m_calculated_mother_ip           = motherParticle.GetDistanceFromVertex( vertex );
  m_calculated_mother_ipchi2       = motherParticle.GetDeviationFromVertex( vertex );
  m_calculated_mother_x            = motherParticle.GetX();
  m_calculated_mother_y            = motherParticle.GetY();
  m_calculated_mother_z            = motherParticle.GetZ();
  m_calculated_mother_px           = motherParticle.GetPx();
  m_calculated_mother_py           = motherParticle.GetPy();
  m_calculated_mother_pz           = motherParticle.GetPz();
  m_calculated_mother_pe           = motherParticle.GetE();
  m_calculated_mother_p            = motherParticle.GetP();
  m_calculated_mother_pt           = motherParticle.GetPt();
  m_calculated_mother_s            = motherParticle.GetS();
  m_calculated_mother_q            = (Int_t) motherParticle.Q();
  m_calculated_mother_eta          = motherParticle.GetEta();
  m_calculated_mother_rapidity     = motherParticle.GetRapidity();
  m_calculated_mother_theta        = motherParticle.GetTheta();
  m_calculated_mother_phi          = motherParticle.GetPhi();
  m_calculated_mother_chi2         = motherParticle.Chi2();
  m_calculated_mother_ndof         = motherParticle.NDF();
  m_calculated_mother_cov          = &motherParticle.CovarianceMatrix()[0];

  KFParticle daughterArray[] = {daughter_1, daughter_2, daughter_3, daughter_4};

  for (int i = 0; i < nTracks; ++i)
  {
      m_calculated_daughter_mass[i]         = daughterArray[i].GetMass();
      m_calculated_daughter_ip[i]           = daughterArray[i].GetDistanceFromVertex( vertex );
      m_calculated_daughter_ipchi2[i]       = daughterArray[i].GetDeviationFromVertex( vertex );
      m_calculated_daughter_x[i]            = daughterArray[i].GetX();
      m_calculated_daughter_y[i]            = daughterArray[i].GetY();
      m_calculated_daughter_z[i]            = daughterArray[i].GetZ();
      m_calculated_daughter_px[i]           = daughterArray[i].GetPx();
      m_calculated_daughter_py[i]           = daughterArray[i].GetPy();
      m_calculated_daughter_pz[i]           = daughterArray[i].GetPz();
      m_calculated_daughter_pe[i]           = daughterArray[i].GetE();
      m_calculated_daughter_p[i]            = daughterArray[i].GetP();
      m_calculated_daughter_pt[i]           = daughterArray[i].GetPt();
      m_calculated_daughter_s[i]            = daughterArray[i].GetS();
      m_calculated_daughter_q[i]            = (Int_t) daughterArray[i].Q();
      m_calculated_daughter_eta[i]          = daughterArray[i].GetEta();
      m_calculated_daughter_rapidity[i]     = daughterArray[i].GetRapidity();
      m_calculated_daughter_theta[i]        = daughterArray[i].GetTheta();
      m_calculated_daughter_phi[i]          = daughterArray[i].GetPhi();
      m_calculated_daughter_chi2[i]         = daughterArray[i].Chi2();
      m_calculated_daughter_ndof[i]         = daughterArray[i].NDF();
      m_calculated_daughter_cov[i]          = &daughterArray[i].CovarianceMatrix()[0];
  }

  m_calculated_vertex_x            = vertex.GetX();
  m_calculated_vertex_y            = vertex.GetY();
  m_calculated_vertex_z            = vertex.GetZ();
  m_calculated_vertex_cov          = &vertex.CovarianceMatrix()[0];

  m_nPVs         = nPVs;
  m_multiplicity = multiplicity;

  m_tree->Fill();
}
