#include <g4main/PHG4Particle.h>
#include <g4eval/SvtxTrackEval.h>
#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxEvalStack.h>
//#include <Acts/Surfaces/Surface.hpp>
//#include "PHActsSourceLinks.h"
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrDefs.h>
#include <mvtx/MvtxDefs.h>
#include <intt/InttDefs.h>
#include <tpc/TpcDefs.h>
#include <KFParticle_nTuple.h>
#include <TTree.h>
#include "KFParticle.h"
#include "KFPVertex.h"
#include <KFParticle_Tools.h>

std::map<std::string, int> Use = 
{
  { "MVTX",  1 },
  { "INTT",  1 },
  { "TPC",   1 },
  { "EMCAL", 0 },
  { "OHCAL", 0 },
  { "IHCAL", 0 }
};

KFParticle_Tools kfpTupleTools;
//For truth matching
SvtxTrackMap *dst_trackmap;
SvtxTrackEval *trackeval;
SvtxClusterEval *clustereval;
SvtxTrack *track;
PHG4Particle* g4particle;
//For detector info

KFParticle_nTuple::KFParticle_nTuple():
    m_truth_matching(false),
    m_detector_info(false),
    m_use_intermediate_name(true),
    m_svtx_evalstack(nullptr)
{} //Constructor

KFParticle_nTuple::~KFParticle_nTuple(){} //Destructor

void KFParticle_nTuple::initializeVariables()
{
  //m_calculated_mother_cov = -99;
}

void KFParticle_nTuple::initializeBranches()
{

  m_tree = new TTree("DecayTree", "DecayTree");
  
  std::string mother_name;
  if (m_mother_name.empty()) mother_name = "mother";
  else mother_name = m_mother_name;

  m_tree->Branch( TString(mother_name) + "_mass",           &m_calculated_mother_mass,            TString(mother_name) + "_mass/F" );
  m_tree->Branch( TString(mother_name) + "_massErr",        &m_calculated_mother_mass_err,        TString(mother_name) + "_massErr/F" );
  m_tree->Branch( TString(mother_name) + "_decayTime",      &m_calculated_mother_decaytime,       TString(mother_name) + "_decayTime/F" );
  m_tree->Branch( TString(mother_name) + "_decayTimeErr",   &m_calculated_mother_decaytime_err,   TString(mother_name) + "_decayTimeErr/F" );
  m_tree->Branch( TString(mother_name) + "_decayLength",    &m_calculated_mother_decaylength,     TString(mother_name) + "_decayLength/F" );
  m_tree->Branch( TString(mother_name) + "_decayLengthErr", &m_calculated_mother_decaylength_err, TString(mother_name) + "_decayLengthErr/F" );
  m_tree->Branch( TString(mother_name) + "_DIRA",           &m_calculated_mother_dira,            TString(mother_name) + "_DIRA/F" );
  m_tree->Branch( TString(mother_name) + "_FDchi2",         &m_calculated_mother_fdchi2,          TString(mother_name) + "_FDchi2/F" );
  m_tree->Branch( TString(mother_name) + "_IP",             &m_calculated_mother_ip,              TString(mother_name) + "_IP/F" );
  m_tree->Branch( TString(mother_name) + "_IPchi2",         &m_calculated_mother_ipchi2,          TString(mother_name) + "_IPchi2/F" );
  m_tree->Branch( TString(mother_name) + "_x",              &m_calculated_mother_x,               TString(mother_name) + "_x/F" );
  m_tree->Branch( TString(mother_name) + "_y",              &m_calculated_mother_y,               TString(mother_name) + "_y/F" );
  m_tree->Branch( TString(mother_name) + "_z",              &m_calculated_mother_z,               TString(mother_name) + "_z/F" );
  m_tree->Branch( TString(mother_name) + "_px",             &m_calculated_mother_px,              TString(mother_name) + "_px/F" );
  m_tree->Branch( TString(mother_name) + "_py",             &m_calculated_mother_py,              TString(mother_name) + "_py/F" );
  m_tree->Branch( TString(mother_name) + "_pz",             &m_calculated_mother_pz,              TString(mother_name) + "_pz/F" );
  m_tree->Branch( TString(mother_name) + "_pE",             &m_calculated_mother_pe,              TString(mother_name) + "_pE/F" );
  m_tree->Branch( TString(mother_name) + "_p",              &m_calculated_mother_p,               TString(mother_name) + "_p/F" );
  m_tree->Branch( TString(mother_name) + "_pErr",           &m_calculated_mother_p_err,           TString(mother_name) + "_pErr/F" );
  m_tree->Branch( TString(mother_name) + "_pT",             &m_calculated_mother_pt,              TString(mother_name) + "_pT/F" );
  m_tree->Branch( TString(mother_name) + "_pTErr",          &m_calculated_mother_pt_err,          TString(mother_name) + "_pTErr/F" );
  m_tree->Branch( TString(mother_name) + "_S",              &m_calculated_mother_s,               TString(mother_name) + "_S/F" );
  m_tree->Branch( TString(mother_name) + "_charge",         &m_calculated_mother_q,               TString(mother_name) + "_charge/I" );
  m_tree->Branch( TString(mother_name) + "_pseudorapidity", &m_calculated_mother_eta,             TString(mother_name) + "_pseudorapidity/F" );
  m_tree->Branch( TString(mother_name) + "_rapidity",       &m_calculated_mother_rapidity,        TString(mother_name) + "_rapidity/F" );
  m_tree->Branch( TString(mother_name) + "_theta",          &m_calculated_mother_theta,           TString(mother_name) + "_theta/F" );
  m_tree->Branch( TString(mother_name) + "_phi",            &m_calculated_mother_phi,             TString(mother_name) + "_phi/F" );
  m_tree->Branch( TString(mother_name) + "_chi2",           &m_calculated_mother_chi2,            TString(mother_name) + "_chi2/F" );
  m_tree->Branch( TString(mother_name) + "_nDoF",           &m_calculated_mother_ndof,            TString(mother_name) + "_nDoF/I" );
  //m_tree->Branch( TString(mother_name) + "_Covariance",      m_calculated_mother_cov,          TString(mother_name) + "_Covariance[21]/F", 21 );

  if ( m_has_intermediates_nTuple )
  {
    for (int i = 0; i < m_num_intermediate_states_nTuple; ++i)
    {
      std::string intermediate_name;
      if (!m_use_intermediate_name) intermediate_name = "intermediate_" + std::to_string(i + 1);
      else intermediate_name = m_intermediate_name_ntuple[i];

      //Note, TBranch will not allow the leaf to contain a forward slash as it is used to define the branch type. Causes problems with J/psi
      std::string fwd_slsh = "/", undrscr = "_"; size_t pos;
      while ((pos = intermediate_name.find(fwd_slsh)) != std::string::npos) intermediate_name.replace(pos, 1, undrscr);
     
      m_tree->Branch( TString(intermediate_name) + "_mass",           &m_calculated_intermediate_mass[i],            TString(intermediate_name) + "_mass/F" );
      m_tree->Branch( TString(intermediate_name) + "_massErr",        &m_calculated_intermediate_mass_err[i],        TString(intermediate_name) + "_massErr/F" );
      m_tree->Branch( TString(intermediate_name) + "_decayTime",      &m_calculated_intermediate_decaytime[i],       TString(intermediate_name) + "_decayTime/F" );
      m_tree->Branch( TString(intermediate_name) + "_decayTimeErr",   &m_calculated_intermediate_decaytime_err[i],   TString(intermediate_name) + "_decayTimeErr/F" );
      m_tree->Branch( TString(intermediate_name) + "_decayLength",    &m_calculated_intermediate_decaylength[i],     TString(intermediate_name) + "_decayLength/F" );
      m_tree->Branch( TString(intermediate_name) + "_decayLengthErr", &m_calculated_intermediate_decaylength_err[i], TString(intermediate_name) + "_decayLengthErr/F" );
      m_tree->Branch( TString(intermediate_name) + "_IP",             &m_calculated_intermediate_ip[i],              TString(intermediate_name) + "_IP/F" );
      m_tree->Branch( TString(intermediate_name) + "_IPchi2",         &m_calculated_intermediate_ipchi2[i],          TString(intermediate_name) + "_IPchi2/F" );
      m_tree->Branch( TString(intermediate_name) + "_x",              &m_calculated_intermediate_x[i],               TString(intermediate_name) + "_x/F" );
      m_tree->Branch( TString(intermediate_name) + "_y",              &m_calculated_intermediate_y[i],               TString(intermediate_name) + "_y/F" );
      m_tree->Branch( TString(intermediate_name) + "_z",              &m_calculated_intermediate_z[i],               TString(intermediate_name) + "_z/F" );
      m_tree->Branch( TString(intermediate_name) + "_px",             &m_calculated_intermediate_px[i],              TString(intermediate_name) + "_px/F" );
      m_tree->Branch( TString(intermediate_name) + "_py",             &m_calculated_intermediate_py[i],              TString(intermediate_name) + "_py/F" );
      m_tree->Branch( TString(intermediate_name) + "_pz",             &m_calculated_intermediate_pz[i],              TString(intermediate_name) + "_pz/F" );
      m_tree->Branch( TString(intermediate_name) + "_pE",             &m_calculated_intermediate_pe[i],              TString(intermediate_name) + "_pE/F" );
      m_tree->Branch( TString(intermediate_name) + "_p",              &m_calculated_intermediate_p[i],               TString(intermediate_name) + "_p/F" );
      m_tree->Branch( TString(intermediate_name) + "_pErr",           &m_calculated_intermediate_p_err[i],           TString(intermediate_name) + "_pErr/F" );
      m_tree->Branch( TString(intermediate_name) + "_pT",             &m_calculated_intermediate_pt[i],              TString(intermediate_name) + "_pT/F" );
      m_tree->Branch( TString(intermediate_name) + "_pTErr",          &m_calculated_intermediate_pt_err[i],          TString(intermediate_name) + "_pTErr/F" );
      m_tree->Branch( TString(intermediate_name) + "_S",              &m_calculated_intermediate_s[i],               TString(intermediate_name) + "_S/F" );
      m_tree->Branch( TString(intermediate_name) + "_charge",         &m_calculated_intermediate_q[i],               TString(intermediate_name) + "_charge/I" );
      m_tree->Branch( TString(intermediate_name) + "_pseudorapidity", &m_calculated_intermediate_eta[i],             TString(intermediate_name) + "_pseudorapidity/F" );
      m_tree->Branch( TString(intermediate_name) + "_rapidity",       &m_calculated_intermediate_rapidity[i],        TString(intermediate_name) + "_rapidity/F" );
      m_tree->Branch( TString(intermediate_name) + "_theta",          &m_calculated_intermediate_theta[i],           TString(intermediate_name) + "_theta/F" );
      m_tree->Branch( TString(intermediate_name) + "_phi",            &m_calculated_intermediate_phi[i],             TString(intermediate_name) + "_phi/F" );
      m_tree->Branch( TString(intermediate_name) + "_chi2",           &m_calculated_intermediate_chi2[i],            TString(intermediate_name) + "_chi2/F" );
      m_tree->Branch( TString(intermediate_name) + "_nDoF",           &m_calculated_intermediate_ndof[i],            TString(intermediate_name) + "_nDoF/I" );
      //m_tree->Branch( TString(intermediate_name) + "_Covariance",      m_calculated_intermediate_cov[i],      TString(intermediate_name) + "_Covariance[21]/F", 21 );     
    }
  }

 for (int i = 0; i < m_num_tracks_nTuple; ++i)
 {
    std::string daughter_number = "track_" + std::to_string(i + 1);

    m_tree->Branch( TString(daughter_number) + "_mass",           &m_calculated_daughter_mass[i],     TString(daughter_number) + "_mass/F" );
    m_tree->Branch( TString(daughter_number) + "_IP",             &m_calculated_daughter_ip[i],       TString(daughter_number) + "_IP/F" );
    m_tree->Branch( TString(daughter_number) + "_IPchi2",         &m_calculated_daughter_ipchi2[i],   TString(daughter_number) + "_IPchi2/F" );
    m_tree->Branch( TString(daughter_number) + "_x",              &m_calculated_daughter_x[i],        TString(daughter_number) + "_x/F" );
    m_tree->Branch( TString(daughter_number) + "_y",              &m_calculated_daughter_y[i],        TString(daughter_number) + "_y/F" );
    m_tree->Branch( TString(daughter_number) + "_z",              &m_calculated_daughter_z[i],        TString(daughter_number) + "_z/F" );
    m_tree->Branch( TString(daughter_number) + "_px",             &m_calculated_daughter_px[i],       TString(daughter_number) + "_px/F" );
    m_tree->Branch( TString(daughter_number) + "_py",             &m_calculated_daughter_py[i],       TString(daughter_number) + "_py/F" );
    m_tree->Branch( TString(daughter_number) + "_pz",             &m_calculated_daughter_pz[i],       TString(daughter_number) + "_pz/F" );
    m_tree->Branch( TString(daughter_number) + "_pE",             &m_calculated_daughter_pe[i],       TString(daughter_number) + "_pE/F" );
    m_tree->Branch( TString(daughter_number) + "_p",              &m_calculated_daughter_p[i],        TString(daughter_number) + "_p/F" );
    m_tree->Branch( TString(daughter_number) + "_pT",             &m_calculated_daughter_pt[i],       TString(daughter_number) + "_pT/F" );
    m_tree->Branch( TString(daughter_number) + "_S",              &m_calculated_daughter_s[i],        TString(daughter_number) + "_S/F" );
    m_tree->Branch( TString(daughter_number) + "_charge",         &m_calculated_daughter_q[i],        TString(daughter_number) + "_charge/I" );
    m_tree->Branch( TString(daughter_number) + "_pseudorapidity", &m_calculated_daughter_eta[i],      TString(daughter_number) + "_pseudorapidity/F" );
    m_tree->Branch( TString(daughter_number) + "_rapidity",       &m_calculated_daughter_rapidity[i], TString(daughter_number) + "_rapidity/F" );
    m_tree->Branch( TString(daughter_number) + "_theta",          &m_calculated_daughter_theta[i],    TString(daughter_number) + "_theta/F" );
    m_tree->Branch( TString(daughter_number) + "_phi",            &m_calculated_daughter_phi[i],      TString(daughter_number) + "_phi/F" );
    m_tree->Branch( TString(daughter_number) + "_chi2",           &m_calculated_daughter_chi2[i],     TString(daughter_number) + "_chi2/F" );
    m_tree->Branch( TString(daughter_number) + "_nDoF",           &m_calculated_daughter_ndof[i],     TString(daughter_number) + "_nDoF/I" );
    m_tree->Branch( TString(daughter_number) + "_trackID",        &m_calculated_daughter_trid[i],     TString(daughter_number) + "_trackID/I" );
    //m_tree->Branch( TString(daughter_number) + "_Covariance",      m_calculated_daughter_cov[i],      TString(daughter_number) + "_Covariance[21]/F", 21 );

    if ( m_truth_matching ) initializeTruthBranches( i ); 
    if ( m_detector_info )  initializeDetectorBranches( i );
 }

  int iter = 0;
  for (int i = 0; i < m_num_tracks_nTuple; ++i) 
  {
    for(int j = 0; j < m_num_tracks_nTuple; ++j)
    {
      if( i < j )
      {
        std::string dca_branch_name = "d" + std::to_string(i + 1) + std::to_string(j + 1) + "_DCA_3D";
        std::string dca_leaf_name = dca_branch_name + "/F";
        m_tree->Branch( dca_branch_name.c_str(), &m_daughter_dca[iter], dca_leaf_name.c_str() );
        ++iter;
      }
    }
  }
  
  m_tree->Branch( "primary_vertex_x",           &m_calculated_vertex_x,   "primary_vertex_x/F" );
  m_tree->Branch( "primary_vertex_y",           &m_calculated_vertex_y,   "primary_vertex_y/F" );
  m_tree->Branch( "primary_vertex_z",           &m_calculated_vertex_z,   "primary_vertex_z/F" );
  //m_tree->Branch( "primary_vertex_Covariance",   m_calculated_vertex_cov, "primary_vertex_Covariance[6]/F", 6 );

  m_tree->Branch( "nPrimaryVertices", &m_nPVs,         "nPrimaryVertices/I" );
  m_tree->Branch( "nEventTracks",     &m_multiplicity, "nEventTracks/I" );

}


void KFParticle_nTuple::fillBranch( PHCompositeNode *topNode,
                                    KFParticle motherParticle, 
                                    KFParticle vertex,
                                    std::vector<KFParticle> daughters,
                                    std::vector<KFParticle> intermediates,
                                    int nPVs, int multiplicity )
{
  KFPVertex *kfpVertex = new KFPVertex; 
  kfpVertex->SetXYZ( vertex.Parameters() );
  kfpVertex->SetCovarianceMatrix( vertex.CovarianceMatrix() ); 
  
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
  m_calculated_mother_p_err        = motherParticle.GetErrP();
  m_calculated_mother_pt           = motherParticle.GetPt();
  m_calculated_mother_pt_err       = motherParticle.GetErrPt();
  m_calculated_mother_s            = motherParticle.GetS();
  m_calculated_mother_q            = (Int_t) motherParticle.Q();
  m_calculated_mother_eta          = motherParticle.GetEta();
  m_calculated_mother_rapidity     = motherParticle.GetRapidity();
  m_calculated_mother_theta        = motherParticle.GetTheta();
  m_calculated_mother_phi          = motherParticle.GetPhi();
  m_calculated_mother_chi2         = motherParticle.GetChi2();
  m_calculated_mother_ndof         = motherParticle.GetNDF();
  //m_calculated_mother_cov          = &motherParticle.CovarianceMatrix()[0];
  
  motherParticle.GetMass( m_calculated_mother_mass, m_calculated_mother_mass_err);
  motherParticle.SetProductionVertex( vertex );
  motherParticle.GetLifeTime( m_calculated_mother_decaytime, m_calculated_mother_decaytime_err );
  motherParticle.GetDecayLength( m_calculated_mother_decaylength, m_calculated_mother_decaylength_err );

  if ( m_has_intermediates_nTuple)
  {
    KFParticle* intermediateArray = &intermediates[0];

    for (int i = 0; i < m_num_intermediate_states_nTuple; ++i)
    {
      m_calculated_intermediate_ip[i]           = intermediateArray[i].GetDistanceFromVertex( vertex );
      m_calculated_intermediate_ipchi2[i]       = intermediateArray[i].GetDeviationFromVertex( vertex );
      m_calculated_intermediate_x[i]            = intermediateArray[i].GetX();
      m_calculated_intermediate_y[i]            = intermediateArray[i].GetY();
      m_calculated_intermediate_z[i]            = intermediateArray[i].GetZ();
      m_calculated_intermediate_px[i]           = intermediateArray[i].GetPx();
      m_calculated_intermediate_py[i]           = intermediateArray[i].GetPy();
      m_calculated_intermediate_pz[i]           = intermediateArray[i].GetPz();
      m_calculated_intermediate_pe[i]           = intermediateArray[i].GetE();
      m_calculated_intermediate_p[i]            = intermediateArray[i].GetP();
      m_calculated_intermediate_p_err[i]        = intermediateArray[i].GetErrP();
      m_calculated_intermediate_pt[i]           = intermediateArray[i].GetPt();
      m_calculated_intermediate_pt_err[i]       = intermediateArray[i].GetErrPt();
      m_calculated_intermediate_s[i]            = intermediateArray[i].GetS();
      m_calculated_intermediate_q[i]            = (Int_t) intermediateArray[i].Q();
      m_calculated_intermediate_eta[i]          = intermediateArray[i].GetEta();
      m_calculated_intermediate_rapidity[i]     = intermediateArray[i].GetRapidity();
      m_calculated_intermediate_theta[i]        = intermediateArray[i].GetTheta();
      m_calculated_intermediate_phi[i]          = intermediateArray[i].GetPhi();
      m_calculated_intermediate_chi2[i]         = intermediateArray[i].GetChi2();
      m_calculated_intermediate_ndof[i]         = intermediateArray[i].GetNDF();
      //m_calculated_intermediate_cov[i]          = &intermediateArray[i].CovarianceMatrix()[0];
      
      intermediateArray[i].GetMass( m_calculated_intermediate_mass[i], m_calculated_intermediate_mass_err[i]);
      intermediateArray[i].SetProductionVertex( motherParticle );
      intermediateArray[i].GetLifeTime( m_calculated_intermediate_decaytime[i], m_calculated_intermediate_decaytime_err[i] );
      intermediateArray[i].GetDecayLength( m_calculated_intermediate_decaylength[i], m_calculated_intermediate_decaylength_err[i] );
    }
  }

  KFParticle temp;
  KFParticle* daughterArray = &daughters[0];

  for( int i = 0; i < m_num_tracks_nTuple; i++ ) //This section of code should rearrange daughter particles by mass
  {
    for( int j = i + 1; j < m_num_tracks_nTuple; j++ )
    {
      if( daughterArray[i].GetMass() > daughterArray[j].GetMass() )
      {
          temp = daughterArray[i];
	  daughterArray[i] = daughterArray[j];
	  daughterArray[j] = temp;
      }
    }
  }

  for ( int i = 0; i < m_num_tracks_nTuple; ++i )
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
      m_calculated_daughter_chi2[i]         = daughterArray[i].GetChi2();
      m_calculated_daughter_ndof[i]         = daughterArray[i].GetNDF();
      m_calculated_daughter_trid[i]         = daughterArray[i].Id(); 
      //m_calculated_daughter_cov[i]          = &daughterArray[i].CovarianceMatrix()[0];

      if ( m_truth_matching ) fillTruthBranch( topNode, daughterArray[i], i ); 
      if ( m_detector_info )  fillDetectorBranch(topNode, daughterArray[i], i );
  }

  int iter = 0;
  for (int i = 0; i < m_num_tracks_nTuple; ++i) 
    for(int j = 0; j < m_num_tracks_nTuple; ++j)
      if( i < j ) { m_daughter_dca[iter] = daughterArray[i].GetDistanceFromParticle( daughterArray[j] ); ++iter; }

  m_calculated_vertex_x            = vertex.GetX();
  m_calculated_vertex_y            = vertex.GetY();
  m_calculated_vertex_z            = vertex.GetZ();
  //m_calculated_vertex_cov          = &vertex.CovarianceMatrix()[0];

  m_nPVs         = nPVs;
  m_multiplicity = multiplicity;

  m_tree->Fill();
}


SvtxTrack* KFParticle_nTuple::getTrack( unsigned int track_id, SvtxTrackMap *trackmap )
{
    SvtxTrack* matched_track = NULL;
    for ( SvtxTrackMap::Iter iter = trackmap->begin(); iter != trackmap->end(); ++iter )
    {
        if (iter->first == track_id) matched_track = iter->second; 
    }
    return matched_track;
}

void KFParticle_nTuple::initializeTruthBranches( int daughter_id )
{
  std::string daughter_number = "daughter_" + std::to_string(daughter_id + 1);

  m_tree->Branch( TString(daughter_number) + "_true_px", &m_true_daughter_px[daughter_id], TString(daughter_number) + "_true_px/F" );
  m_tree->Branch( TString(daughter_number) + "_true_py", &m_true_daughter_py[daughter_id], TString(daughter_number) + "_true_py/F" );
  m_tree->Branch( TString(daughter_number) + "_true_pz", &m_true_daughter_pz[daughter_id], TString(daughter_number) + "_true_pz/F" );
}

void  KFParticle_nTuple::fillTruthBranch( PHCompositeNode *topNode, KFParticle daughter, int daughter_id )
{
    if (!m_svtx_evalstack)
    {
      m_svtx_evalstack = new SvtxEvalStack(topNode);
      m_svtx_evalstack->set_strict(true);
      trackeval = m_svtx_evalstack->get_track_eval();
      clustereval = m_svtx_evalstack->get_cluster_eval();
    }

    dst_trackmap = findNode::getClass<SvtxTrackMap>( topNode, "SvtxTrackMap" );
    track = getTrack( daughter.Id(), dst_trackmap );

    SvtxTrack::ConstClusterKeyIter iter = track->begin_cluster_keys();
    TrkrDefs::cluskey clusKey = *iter;

    if (!clusKey) g4particle = trackeval->max_truth_particle_by_nclusters( track );
    else g4particle = clustereval->max_truth_particle_by_energy(clusKey);

    m_true_daughter_px[ daughter_id ] = (Float_t) g4particle->get_px();
    m_true_daughter_py[ daughter_id ] = (Float_t) g4particle->get_py();
    m_true_daughter_pz[ daughter_id ] = (Float_t) g4particle->get_pz();
}

void KFParticle_nTuple::initializeDetectorBranches( int daughter_id )
{
    std::string daughter_number = "daughter_" + std::to_string(daughter_id + 1);

    m_tree->Branch( TString(daughter_number) + "_local_x", &detector_local_x[ daughter_id ] );
    m_tree->Branch( TString(daughter_number) + "_local_y", &detector_local_y[ daughter_id ] );
    m_tree->Branch( TString(daughter_number) + "_local_z", &detector_local_z[ daughter_id ] );
    m_tree->Branch( TString(daughter_number) +"_layer", &detector_layer[ daughter_id ] );

    for ( auto const& subdetector : Use )
    {
      if ( subdetector.second ) initializeSubDetectorBranches( subdetector.first, daughter_id );
    }
}


void KFParticle_nTuple::initializeSubDetectorBranches( std::string detectorName, int daughter_id )
{
    std::string daughter_number = "daughter_" + std::to_string(daughter_id + 1);

    if (detectorName == "MVTX" ) m_tree->Branch( TString(daughter_number) + "_" + TString(detectorName) + "_staveID", &mvtx_staveID[ daughter_id ] );
    if (detectorName == "MVTX" ) m_tree->Branch( TString(daughter_number) + "_" + TString(detectorName) + "_chipID", &mvtx_chipID[ daughter_id ] );
    if (detectorName == "INTT" ) m_tree->Branch( TString(daughter_number) + "_" + TString(detectorName) + "_ladderZID", &intt_ladderZID[ daughter_id ] );
    if (detectorName == "INTT" ) m_tree->Branch( TString(daughter_number) + "_" + TString(detectorName) + "_ladderPhiID", &intt_ladderPhiID[ daughter_id ] );
    if (detectorName == "TPC" )  m_tree->Branch( TString(daughter_number) + "_" + TString(detectorName) + "_sectorID", &tpc_sectorID[ daughter_id ] );
    if (detectorName == "TPC" )  m_tree->Branch( TString(daughter_number) + "_" + TString(detectorName) + "_side", &tpc_side[ daughter_id ] );
}

void KFParticle_nTuple::fillDetectorBranch( PHCompositeNode *topNode, KFParticle daughter, int daughter_id )
{
    dst_clustermap = findNode::getClass<TrkrClusterContainer>( topNode, "TRKR_CLUSTER");
    dst_trackmap = findNode::getClass<SvtxTrackMap>( topNode, "SvtxTrackMap" );
    track = getTrack( daughter.Id(), dst_trackmap );

    for (SvtxTrack::ConstClusterKeyIter iter = track->begin_cluster_keys(); iter != track->end_cluster_keys(); ++iter)
    {
       TrkrDefs::cluskey clusKey = *iter;
       TrkrCluster *cluster = dst_clustermap->findCluster(clusKey);
       const unsigned int trkrId = TrkrDefs::getTrkrId(clusKey);

       detector_local_x[ daughter_id ].push_back( cluster->getX() );
       detector_local_y[ daughter_id ].push_back( cluster->getY() );
       detector_local_z[ daughter_id ].push_back( cluster->getZ() );
       detector_layer[ daughter_id ].push_back( TrkrDefs::getLayer(clusKey) );
       unsigned int staveId, chipId, ladderZId, ladderPhiId, sectorId, side;     
       staveId = chipId = ladderZId = ladderPhiId = sectorId = side = -99;

      if ( Use["MVTX"] && trkrId == TrkrDefs::mvtxId )
      {
        staveId = MvtxDefs::getStaveId(clusKey);
        chipId = MvtxDefs::getChipId(clusKey);
      } 
      else if (  Use["INTT"] && trkrId == TrkrDefs::inttId )
      {
        ladderZId = InttDefs::getLadderZId(clusKey);
        ladderPhiId = InttDefs::getLadderPhiId(clusKey);
      }
      else if (  Use["TPC"] && trkrId == TrkrDefs::tpcId )
      {
        sectorId = TpcDefs::getSectorId(clusKey);
        side = TpcDefs::getSide(clusKey);
      }
      
      mvtx_staveID[ daughter_id ].push_back( staveId );
      mvtx_chipID[ daughter_id ].push_back( chipId );
      intt_ladderZID[ daughter_id ].push_back( ladderZId );
      intt_ladderPhiID[ daughter_id ].push_back( ladderPhiId );
      tpc_sectorID[ daughter_id ].push_back( sectorId );
      tpc_side[ daughter_id ].push_back( side );

   }

}
