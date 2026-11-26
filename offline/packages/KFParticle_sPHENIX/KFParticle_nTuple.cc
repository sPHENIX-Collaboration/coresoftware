#include "KFParticle_nTuple.h"
#include "KFParticle_util.h"

#include <ffaobjects/EventHeader.h>
#include <ffarawobjects/Gl1Packet.h>
#include <fun4all/Fun4AllHistoManager.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE
#include <qautils/QAHistManagerDef.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackAnalysisUtils.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase/TrkrDefs.h>

#include <KFParticle.h>
#include <KFVertex.h>

#include <centrality/CentralityInfo.h>

#include <Rtypes.h>
#include <TString.h>  // for TString, operator+
#include <TTree.h>

#include <Eigen/Dense>

#include <algorithm>  // for max
#include <cmath>
#include <cstdlib>  // for abs, size_t
#include <map>      // for map, _Rb_tree_iterator, map<>:...
#include <utility>  // for swap

class PHCompositeNode;
class PHNode;

float TempError;

void KFParticle_nTuple::initializeVariables()
{
  // m_calculated_mother_cov = -99;
}

void KFParticle_nTuple::initializeBranches(PHCompositeNode* topNode)
{
  // KFParticle_Tools::init_dEdx_fits(); //Cant do this! Two trees open at once!

  delete m_tree;
  m_tree = new TTree("DecayTree", "DecayTree");
  m_tree->OptimizeBaskets();
  m_tree->SetAutoSave(-5e6);  // Save the output file every 5MB

  // Note, TBranch will not allow the leaf to contain a forward slash as it is used to define the branch type
  const std::map<std::string, std::string> forbiddenStrings = {
      {"/", "_"}, {"(", "_"}, {")", ""}, {"+", "plus"}, {"-", "minus"}, {"*", "star"}};
  auto sanitize = [&forbiddenStrings](std::string name)
  {
    size_t pos;
    for (const auto& [badString, goodString] : forbiddenStrings)
    {
      while ((pos = name.find(badString)) != std::string::npos)
      {
        name.replace(pos, 1, goodString);
      }
    }
    return name;
  };

  std::string mother_name = sanitize(m_config.motherName.empty() ? "mother" : m_config.motherName);
  //setup only the momentum calibration branches
  
  if (m_config.lightweightNTuple)
  {
    initializeMomentumCalibrationBranches(mother_name);
    // Hand the tree to the QA histogram manager, which then owns it (and detaches it from any TFile)
    if (m_config.lightweightNTupleToQA && !QAHistManagerDef::getHistoManager()->registerHisto(m_qa_tree_path, m_tree))
    {
      std::cout << "KFParticle_nTuple: could not register " << m_qa_tree_path << " with the QA histogram manager" << std::endl;
    }
    return;
  }

  m_tree->Branch(TString(mother_name) + "_mass", &m_calculated_mother_mass, TString(mother_name) + "_mass/F");
  m_tree->Branch(TString(mother_name) + "_massErr", &m_calculated_mother_mass_err, TString(mother_name) + "_massErr/F");
  if (m_config.constrainToPrimaryVertex)
  {
    m_tree->Branch(TString(mother_name) + "_decayTime", &m_calculated_mother_decaytime, TString(mother_name) + "_decayTime/F");
    m_tree->Branch(TString(mother_name) + "_decayTimeErr", &m_calculated_mother_decaytime_err, TString(mother_name) + "_decayTimeErr/F");
    m_tree->Branch(TString(mother_name) + "_decayLength", &m_calculated_mother_decaylength, TString(mother_name) + "_decayLength/F");
    m_tree->Branch(TString(mother_name) + "_decayLengthErr", &m_calculated_mother_decaylength_err, TString(mother_name) + "_decayLengthErr/F");
    m_tree->Branch(TString(mother_name) + "_decayLength_xy", &m_calculated_mother_decaylength_xy, TString(mother_name) + "_decayLength_xy/F");
    m_tree->Branch(TString(mother_name) + "_decayLengthErr_xy", &m_calculated_mother_decaylength_xy_err, TString(mother_name) + "_decayLengthErr_xy/F");
    m_tree->Branch(TString(mother_name) + "_DIRA", &m_calculated_mother_dira, TString(mother_name) + "_DIRA/F");
    m_tree->Branch(TString(mother_name) + "_DIRA_xy", &m_calculated_mother_dira_xy, TString(mother_name) + "_DIRA_xy/F");
    m_tree->Branch(TString(mother_name) + "_FDchi2", &m_calculated_mother_fdchi2, TString(mother_name) + "_FDchi2/F");
    m_tree->Branch(TString(mother_name) + "_PV_DCA", &m_calculated_mother_PV_dca, TString(mother_name) + "_PV_DCA/F");
    m_tree->Branch(TString(mother_name) + "_PV_DCA_StdDev", &m_calculated_mother_PV_dca_sig, TString(mother_name) + "_PV_DCA_StdDev/F");
    m_tree->Branch(TString(mother_name) + "_PV_DCA_Err", &m_calculated_mother_PV_dca_err, TString(mother_name) + "_PV_DCA_Err/F");
    m_tree->Branch(TString(mother_name) + "_PV_DCA_xy", &m_calculated_mother_PV_dca_xy, TString(mother_name) + "_PV_DCA_xy/F");
  }
  if (m_config.allPVInfo)
  {
    m_tree->Branch(TString(mother_name) + "_PV_DCA_allPV", &allPV_mother_PV_DCA);
    m_tree->Branch(TString(mother_name) + "_PV_DCA_StdDev_allPV", &allPV_mother_PV_DCA_StdDev);
  }
  m_tree->Branch(TString(mother_name) + "_x", &m_calculated_mother_x, TString(mother_name) + "_x/F");
  m_tree->Branch(TString(mother_name) + "_y", &m_calculated_mother_y, TString(mother_name) + "_y/F");
  m_tree->Branch(TString(mother_name) + "_z", &m_calculated_mother_z, TString(mother_name) + "_z/F");
  m_tree->Branch(TString(mother_name) + "_px", &m_calculated_mother_px, TString(mother_name) + "_px/F");
  m_tree->Branch(TString(mother_name) + "_py", &m_calculated_mother_py, TString(mother_name) + "_py/F");
  m_tree->Branch(TString(mother_name) + "_pz", &m_calculated_mother_pz, TString(mother_name) + "_pz/F");
  m_tree->Branch(TString(mother_name) + "_pE", &m_calculated_mother_pe, TString(mother_name) + "_pE/F");
  m_tree->Branch(TString(mother_name) + "_p", &m_calculated_mother_p, TString(mother_name) + "_p/F");
  m_tree->Branch(TString(mother_name) + "_pErr", &m_calculated_mother_p_err, TString(mother_name) + "_pErr/F");
  m_tree->Branch(TString(mother_name) + "_pT", &m_calculated_mother_pt, TString(mother_name) + "_pT/F");
  m_tree->Branch(TString(mother_name) + "_pTErr", &m_calculated_mother_pt_err, TString(mother_name) + "_pTErr/F");
  m_tree->Branch(TString(mother_name) + "_charge", &m_calculated_mother_q, TString(mother_name) + "_charge/B");
  m_tree->Branch(TString(mother_name) + "_pseudorapidity", &m_calculated_mother_eta, TString(mother_name) + "_pseudorapidity/F");
  m_tree->Branch(TString(mother_name) + "_rapidity", &m_calculated_mother_rapidity, TString(mother_name) + "_rapidity/F");
  m_tree->Branch(TString(mother_name) + "_theta", &m_calculated_mother_theta, TString(mother_name) + "_theta/F");
  m_tree->Branch(TString(mother_name) + "_phi", &m_calculated_mother_phi, TString(mother_name) + "_phi/F");
  m_tree->Branch(TString(mother_name) + "_vertex_volume", &m_calculated_mother_v, TString(mother_name) + "_vertex_volume/F");
  m_tree->Branch(TString(mother_name) + "_chi2", &m_calculated_mother_chi2, TString(mother_name) + "_chi2/F");
  m_tree->Branch(TString(mother_name) + "_nDoF", &m_calculated_mother_ndof, TString(mother_name) + "_nDoF/i");
  m_tree->Branch(TString(mother_name) + "_PDG_ID", &m_calculated_mother_pdgID, TString(mother_name) + "_PDG_ID/I");
  m_tree->Branch(TString(mother_name) + "_Covariance", &m_calculated_mother_cov, TString(mother_name) + "_Covariance[21]/F", 21);

  std::vector<std::string> intermediateNameMapping;  // What intermediate is associate to what track
  if (m_config.hasIntermediates)
  {
    for (int i = 0; i < m_config.numIntermediateStates; ++i)
    {
      std::string intermediate_name = sanitize(m_config.intermediateNames[i]);

      m_tree->Branch(TString(intermediate_name) + "_mass", &m_calculated_intermediate_mass[i], TString(intermediate_name) + "_mass/F");
      m_tree->Branch(TString(intermediate_name) + "_massErr", &m_calculated_intermediate_mass_err[i], TString(intermediate_name) + "_massErr/F");
      m_tree->Branch(TString(intermediate_name) + "_decayTime", &m_calculated_intermediate_decaytime[i], TString(intermediate_name) + "_decayTime/F");
      m_tree->Branch(TString(intermediate_name) + "_decayTimeErr", &m_calculated_intermediate_decaytime_err[i], TString(intermediate_name) + "_decayTimeErr/F");
      m_tree->Branch(TString(intermediate_name) + "_decayLength", &m_calculated_intermediate_decaylength[i], TString(intermediate_name) + "_decayLength/F");
      m_tree->Branch(TString(intermediate_name) + "_decayLengthErr", &m_calculated_intermediate_decaylength_err[i], TString(intermediate_name) + "_decayLengthErr/F");
      m_tree->Branch(TString(intermediate_name) + "_decayLength_xy", &m_calculated_intermediate_decaylength_xy[i], TString(intermediate_name) + "_decayLength_xy/F");
      m_tree->Branch(TString(intermediate_name) + "_decayLengthErr_xy", &m_calculated_intermediate_decaylength_xy_err[i], TString(intermediate_name) + "_decayLengthErr_xy/F");
      m_tree->Branch(TString(intermediate_name) + "_DIRA", &m_calculated_intermediate_dira[i], TString(intermediate_name) + "_DIRA/F");
      m_tree->Branch(TString(intermediate_name) + "_FDchi2", &m_calculated_intermediate_fdchi2[i], TString(intermediate_name) + "_FDchi2/F");
      if (m_config.constrainToPrimaryVertex)
      {
        m_tree->Branch(TString(intermediate_name) + "_PV_DCA", &m_calculated_intermediate_PV_dca[i], TString(intermediate_name) + "_PV_DCA/F");
        m_tree->Branch(TString(intermediate_name) + "_PV_DCA_StdDev", &m_calculated_intermediate_PV_dca_sig[i], TString(intermediate_name) + "_PV_DCA_StdDev/F");
        m_tree->Branch(TString(intermediate_name) + "_PV_DCA_Err", &m_calculated_intermediate_PV_dca_err[i], TString(intermediate_name) + "_PV_DCA_Err/F");
        m_tree->Branch(TString(intermediate_name) + "_PV_DCA_xy", &m_calculated_intermediate_PV_dca_xy[i], TString(intermediate_name) + "_PV_DCA_xy/F");
      }
      if (m_config.allPVInfo)
      {
        m_tree->Branch(TString(intermediate_name) + "_PV_DCA_allPV", &allPV_intermediates_PV_DCA[i]);
        m_tree->Branch(TString(intermediate_name) + "_PV_DCA_StdDev_allPV", &allPV_intermediates_PV_DCA_StdDev[i]);
      }
      m_tree->Branch(TString(intermediate_name) + "_x", &m_calculated_intermediate_x[i], TString(intermediate_name) + "_x/F");
      m_tree->Branch(TString(intermediate_name) + "_y", &m_calculated_intermediate_y[i], TString(intermediate_name) + "_y/F");
      m_tree->Branch(TString(intermediate_name) + "_z", &m_calculated_intermediate_z[i], TString(intermediate_name) + "_z/F");
      m_tree->Branch(TString(intermediate_name) + "_px", &m_calculated_intermediate_px[i], TString(intermediate_name) + "_px/F");
      m_tree->Branch(TString(intermediate_name) + "_py", &m_calculated_intermediate_py[i], TString(intermediate_name) + "_py/F");
      m_tree->Branch(TString(intermediate_name) + "_pz", &m_calculated_intermediate_pz[i], TString(intermediate_name) + "_pz/F");
      m_tree->Branch(TString(intermediate_name) + "_pE", &m_calculated_intermediate_pe[i], TString(intermediate_name) + "_pE/F");
      m_tree->Branch(TString(intermediate_name) + "_p", &m_calculated_intermediate_p[i], TString(intermediate_name) + "_p/F");
      m_tree->Branch(TString(intermediate_name) + "_pErr", &m_calculated_intermediate_p_err[i], TString(intermediate_name) + "_pErr/F");
      m_tree->Branch(TString(intermediate_name) + "_pT", &m_calculated_intermediate_pt[i], TString(intermediate_name) + "_pT/F");
      m_tree->Branch(TString(intermediate_name) + "_pTErr", &m_calculated_intermediate_pt_err[i], TString(intermediate_name) + "_pTErr/F");
      m_tree->Branch(TString(intermediate_name) + "_charge", &m_calculated_intermediate_q[i], TString(intermediate_name) + "_charge/B");
      m_tree->Branch(TString(intermediate_name) + "_pseudorapidity", &m_calculated_intermediate_eta[i], TString(intermediate_name) + "_pseudorapidity/F");
      m_tree->Branch(TString(intermediate_name) + "_rapidity", &m_calculated_intermediate_rapidity[i], TString(intermediate_name) + "_rapidity/F");
      m_tree->Branch(TString(intermediate_name) + "_theta", &m_calculated_intermediate_theta[i], TString(intermediate_name) + "_theta/F");
      m_tree->Branch(TString(intermediate_name) + "_phi", &m_calculated_intermediate_phi[i], TString(intermediate_name) + "_phi/F");
      m_tree->Branch(TString(intermediate_name) + "_vertex_volume", &m_calculated_intermediate_v[i], TString(intermediate_name) + "_vertex_volume/F");
      m_tree->Branch(TString(intermediate_name) + "_chi2", &m_calculated_intermediate_chi2[i], TString(intermediate_name) + "_chi2/F");
      m_tree->Branch(TString(intermediate_name) + "_nDoF", &m_calculated_intermediate_ndof[i], TString(intermediate_name) + "_nDoF/i");
      m_tree->Branch(TString(intermediate_name) + "_PDG_ID", &m_calculated_intermediate_pdgID[i], TString(intermediate_name) + "_PDG_ID/I");
      m_tree->Branch(TString(intermediate_name) + "_Covariance", &m_calculated_intermediate_cov[i], TString(intermediate_name) + "_Covariance[21]/F", 21);

      for (int j = 0; j < m_config.numTracksFromIntermediate[i]; ++j)
      {
        intermediateNameMapping.push_back(intermediate_name + "_");
      }
    }
  }

  int num_intermediate_tracks = 0;
  for (int i = 0; i < m_config.numIntermediateStates; ++i)
  {
    num_intermediate_tracks += m_config.numTracksFromIntermediate[i];
  }

  for (int i = 0; i < m_config.numTracks; ++i)
  {
    std::string daughter_number = "track_" + std::to_string(i + 1);

    if (m_config.hasIntermediates && i < num_intermediate_tracks)
    {
      daughter_number.insert(0, intermediateNameMapping[i]);
    }

    m_tree->Branch(TString(daughter_number) + "_mass", &m_calculated_daughter_mass[i], TString(daughter_number) + "_mass/F");
    if (m_config.constrainToPrimaryVertex)
    {
      m_tree->Branch(TString(daughter_number) + "_PV_DCA", &m_calculated_daughter_PV_dca[i], TString(daughter_number) + "_PV_DCA/F");
      m_tree->Branch(TString(daughter_number) + "_PV_DCA_Err", &m_calculated_daughter_PV_dca_err[i], TString(daughter_number) + "_PV_DCA_Err/F");
      m_tree->Branch(TString(daughter_number) + "_PV_DCA_xy", &m_calculated_daughter_PV_dca_xy[i], TString(daughter_number) + "_PV_DCA_xy/F");
      m_tree->Branch(TString(daughter_number) + "_PV_DCA_sig", &m_calculated_daughter_PV_dca_sig[i], TString(daughter_number) + "_PV_DCA_sig/F");
      m_tree->Branch(TString(daughter_number) + "_PV_DCA_sig_xy", &m_calculated_daughter_PV_dca_xy_sig[i], TString(daughter_number) + "_PV_DCA_sig_xy/F");
    }
    if (m_config.allPVInfo)
    {
      m_tree->Branch(TString(daughter_number) + "_PV_DCA_allPV", &allPV_daughter_PV_DCA[i]);
      m_tree->Branch(TString(daughter_number) + "_PV_DCA_StdDev_allPV", &allPV_daughter_PV_DCA_StdDev[i]);
    }
    m_tree->Branch(TString(daughter_number) + "_x", &m_calculated_daughter_x[i], TString(daughter_number) + "_x/F");
    m_tree->Branch(TString(daughter_number) + "_y", &m_calculated_daughter_y[i], TString(daughter_number) + "_y/F");
    m_tree->Branch(TString(daughter_number) + "_z", &m_calculated_daughter_z[i], TString(daughter_number) + "_z/F");
    m_tree->Branch(TString(daughter_number) + "_px", &m_calculated_daughter_px[i], TString(daughter_number) + "_px/F");
    m_tree->Branch(TString(daughter_number) + "_py", &m_calculated_daughter_py[i], TString(daughter_number) + "_py/F");
    m_tree->Branch(TString(daughter_number) + "_pz", &m_calculated_daughter_pz[i], TString(daughter_number) + "_pz/F");
    m_tree->Branch(TString(daughter_number) + "_pE", &m_calculated_daughter_pe[i], TString(daughter_number) + "_pE/F");
    m_tree->Branch(TString(daughter_number) + "_p", &m_calculated_daughter_p[i], TString(daughter_number) + "_p/F");
    m_tree->Branch(TString(daughter_number) + "_pErr", &m_calculated_daughter_p_err[i], TString(daughter_number) + "_pErr/F");
    m_tree->Branch(TString(daughter_number) + "_pT", &m_calculated_daughter_pt[i], TString(daughter_number) + "_pT/F");
    m_tree->Branch(TString(daughter_number) + "_pTErr", &m_calculated_daughter_pt_err[i], TString(daughter_number) + "_pTErr/F");
    m_tree->Branch(TString(daughter_number) + "_jT", &m_calculated_daughter_jt[i], TString(daughter_number) + "_jT/F");
    m_tree->Branch(TString(daughter_number) + "_charge", &m_calculated_daughter_q[i], TString(daughter_number) + "_charge/B");
    m_tree->Branch(TString(daughter_number) + "_bunch_crossing", &m_calculated_daughter_bunch_crossing[i], TString(daughter_number) + "_bunch_crossing/I");
    m_tree->Branch(TString(daughter_number) + "_pseudorapidity", &m_calculated_daughter_eta[i], TString(daughter_number) + "_pseudorapidity/F");
    m_tree->Branch(TString(daughter_number) + "_rapidity", &m_calculated_daughter_rapidity[i], TString(daughter_number) + "_rapidity/F");
    m_tree->Branch(TString(daughter_number) + "_theta", &m_calculated_daughter_theta[i], TString(daughter_number) + "_theta/F");
    m_tree->Branch(TString(daughter_number) + "_phi", &m_calculated_daughter_phi[i], TString(daughter_number) + "_phi/F");
    m_tree->Branch(TString(daughter_number) + "_chi2", &m_calculated_daughter_chi2[i], TString(daughter_number) + "_chi2/F");
    m_tree->Branch(TString(daughter_number) + "_nDoF", &m_calculated_daughter_ndof[i], TString(daughter_number) + "_nDoF/i");
    m_tree->Branch(TString(daughter_number) + "_track_ID", &m_calculated_daughter_trid[i], TString(daughter_number) + "_track_ID/I");
    m_tree->Branch(TString(daughter_number) + "_PDG_ID", &m_calculated_daughter_pdgID[i], TString(daughter_number) + "_PDG_ID/I");
    m_tree->Branch(TString(daughter_number) + "_Covariance", &m_calculated_daughter_cov[i], TString(daughter_number) + "_Covariance[21]/F", 21);
    if (m_config.dEdxInfo)
    {
      m_tree->Branch(TString(daughter_number) + "_calculated_dEdx", &m_calculated_daughter_dedx[i], TString(daughter_number) + "_calculated_dEdx/F");  // m_config.dEdxInfo defaults to false; run get_dEdx_info() to change this
    }
    // m_tree->Branch(TString(daughter_number) + "_expected_pion_dEdx", &m_calculated_daughter_expected_dedx_pion[i], TString(daughter_number) + "_expected_pion_dEdx/F");
    // m_tree->Branch(TString(daughter_number) + "_expected_kaon_dEdx", &m_calculated_daughter_expected_dedx_kaon[i], TString(daughter_number) + "_expected_kaon_dEdx/F");
    // m_tree->Branch(TString(daughter_number) + "_expected_proton_dEdx", &m_calculated_daughter_expected_dedx_proton[i], TString(daughter_number) + "_expected_proton_dEdx/F");
    if (m_config.caloInfo)
    {
      m_calo.initializeBranches(m_tree, i, daughter_number);
    }
    if (m_config.truthMatching)
    {
      m_truth.initializeBranches(m_tree, i, daughter_number, m_config.constrainToPrimaryVertex);
    }
    if (m_config.detectorInfo)
    {
      m_detector.initializeBranches(m_tree, i, daughter_number);
    }
  }

  int iter = 0;
  for (int i = 0; i < m_config.numTracks; ++i)
  {
    for (int j = i + 1; j < m_config.numTracks; ++j)
    {
      const std::string branch_name_base = "track_" + std::to_string(i + 1) + "_track_" + std::to_string(j + 1);
      std::string dca_branch_name = branch_name_base + "_DCA";
      std::string dca_leaf_name = dca_branch_name + "/F";
      m_tree->Branch(dca_branch_name.c_str(), &m_daughter_dca[iter], dca_leaf_name.c_str());

      std::string dca_branch_name_xy = dca_branch_name + "_xy";
      std::string dca_leaf_name_xy = dca_branch_name_xy + "/F";
      m_tree->Branch(dca_branch_name_xy.c_str(), &m_daughter_dca_xy[iter], dca_leaf_name_xy.c_str());

      std::string dca_sig_branch_name = branch_name_base + "_DCA_sig";
      std::string dca_sig_leaf_name = dca_sig_branch_name + "/F";
      m_tree->Branch(dca_sig_branch_name.c_str(), &m_daughter_dca_sig[iter], dca_sig_leaf_name.c_str());

      std::string dca_sig_branch_name_xy = dca_sig_branch_name + "_xy";
      std::string dca_sig_leaf_name_xy = dca_sig_branch_name_xy + "/F";
      m_tree->Branch(dca_sig_branch_name_xy.c_str(), &m_daughter_dca_sig_xy[iter], dca_sig_leaf_name_xy.c_str());

      ++iter;
    }
  }

  if (m_config.constrainToPrimaryVertex)
  {
    m_tree->Branch("primary_vertex_x", &m_calculated_vertex_x, "primary_vertex_x/F");
    m_tree->Branch("primary_vertex_y", &m_calculated_vertex_y, "primary_vertex_y/F");
    m_tree->Branch("primary_vertex_z", &m_calculated_vertex_z, "primary_vertex_z/F");
    m_tree->Branch("primary_vertex_volume", &m_calculated_vertex_v, "primary_vertex_volume/F");
    m_tree->Branch("primary_vertex_chi2", &m_calculated_vertex_chi2, "primary_vertex_chi2/F");
    m_tree->Branch("primary_vertex_nDoF", &m_calculated_vertex_ndof, "primary_vertex_nDoF/i");
    m_tree->Branch("primary_vertex_ID", &m_calculated_vertex_ID, "primary_vertex_ID/I");
    // m_tree->Branch( "primary_vertex_Covariance",   m_calculated_vertex_cov, "primary_vertex_Covariance[6]/F", 6 );
    m_tree->Branch("primary_vertex_Covariance", &m_calculated_vertex_cov, "primary_vertex_Covariance[6]/F", 6);
  }

  if (m_config.useCentrality)
  {
    m_tree->Branch("centrality_MBD", &m_centrality_mbd);
  }

  if (m_config.allPVInfo)
  {
    m_tree->Branch("all_primary_vertex_x", &allPV_x);
    m_tree->Branch("all_primary_vertex_y", &allPV_y);
    m_tree->Branch("all_primary_vertex_z", &allPV_z);
  }

  m_tree->Branch("secondary_vertex_mass_pionPID", &m_sv_mass, "secondary_vertex_mass_pionPID/F");

  m_tree->Branch("nPrimaryVerticesOfBC", &m_nPVs, "nPrimaryVerticesOfBC/I");
  m_tree->Branch("nTracksOfBC", &m_multiplicity, "nTracksOfBC/I");
  m_tree->Branch("nSiSeedMultiplicity", &m_ncharged_siseed_multiplicity, "nSiSeedMultiplicity/I");
  m_tree->Branch("nTracksOfVertex", &m_nTracksOfVertex, "nTracksOfVertex/I");

  m_tree->Branch("runNumber", &m_runNumber, "runNumber/I");
  m_tree->Branch("eventNumber", &m_evtNumber, "eventNumber/I");
  m_tree->Branch("Collision_BCO", &m_bco, "Collision_BCO/L"); //already there, this is shifted BCO
  m_tree->Branch("GL1_BCO", &m_event_bco, "GL1_BCO/L"); //adding for the current event BCO, not shifted
  m_tree->Branch("last_GL1_BCO", &m_last_event_bco, "last_GL1_BCO/L"); //BCO for the last event 

  if (m_config.triggerInfo)
  {
    m_trigger.initializeBranches(topNode, m_tree);
  }
}

void KFParticle_nTuple::fillBranch(PHCompositeNode* topNode,
                                   KFParticle motherParticle,
                                   const KFParticle& vertex_fillbranch,
                                   std::vector<KFParticle> daughters,
                                   std::vector<KFParticle> intermediates)
{
  KFParticle* daughterArray = daughters.data();

  int num_tracks_used_by_intermediates = 0;

  // Order the n tracks starting at num_tracks_used_by_intermediates (by PDG ID, or by mass for charge conjugates)
  auto sortTracks = [&](int n)
  {
    for (int i = 0; i < n; ++i)
    {
      for (int j = i + 1; j < n; ++j)
      {
        int particleAElement = i + num_tracks_used_by_intermediates;
        int particleBElement = j + num_tracks_used_by_intermediates;
        int particleAPID = daughterArray[particleAElement].GetPDG();
        int particleBPID = daughterArray[particleBElement].GetPDG();

        bool switchTrackPosition;
        if (m_config.chargeConjugate)
        {
          float daughterA_mass = KFParticle_util::getParticleMass(particleAPID);
          float daughterB_mass = KFParticle_util::getParticleMass(particleBPID);
          switchTrackPosition = daughterA_mass > daughterB_mass;
        }
        else
        {
          switchTrackPosition = particleAPID > particleBPID;
        }
        if (switchTrackPosition)
        {
          std::swap(daughterArray[particleAElement], daughterArray[particleBElement]);
        }
      }
    }
  };

  for (int k = 0; k < m_config.numIntermediateStates; ++k)  // Rearrange tracks from intermediate states
  {
    sortTracks(m_config.numTracksFromIntermediate[k]);
    num_tracks_used_by_intermediates += m_config.numTracksFromIntermediate[k];
  }

  int num_remaining_tracks = m_config.numTracks - num_tracks_used_by_intermediates;

  sortTracks(num_remaining_tracks);

  if (m_config.extrapolateTracksToSV)
  {
    for (auto &daughter : daughters)
    {
      daughter.SetProductionVertex(motherParticle);
    }
  }
  //This is the setup for filling only the momentum calibration branches
  //It is placed here after the daughters are prepared but before the branches are filled
  if (m_config.lightweightNTuple)
  {
    fillMomentumCalibrationBranch(topNode, motherParticle, daughterArray);
    return;
  }

  if (m_config.constrainToPrimaryVertex)
  {
    m_calculated_mother_dira = KFParticle_util::eventDIRA(motherParticle, vertex_fillbranch);
    m_calculated_mother_dira_xy = KFParticle_util::eventDIRA(motherParticle, vertex_fillbranch, false);
    m_calculated_mother_fdchi2 = KFParticle_util::flightDistanceChi2(motherParticle, vertex_fillbranch);
    m_calculated_mother_PV_dca = motherParticle.GetDistanceFromVertex(vertex_fillbranch);
    m_calculated_mother_PV_dca_sig = motherParticle.GetDeviationFromVertex(vertex_fillbranch);
    m_calculated_mother_PV_dca_err = m_calculated_mother_PV_dca / std::sqrt(m_calculated_mother_PV_dca_sig);
    m_calculated_mother_PV_dca_xy = motherParticle.GetDistanceFromVertexXY(vertex_fillbranch);
  }
  m_calculated_mother_x = motherParticle.GetX();
  m_calculated_mother_y = motherParticle.GetY();
  m_calculated_mother_z = motherParticle.GetZ();
  m_calculated_mother_px = motherParticle.GetPx();
  m_calculated_mother_py = motherParticle.GetPy();
  m_calculated_mother_pz = motherParticle.GetPz();
  m_calculated_mother_pe = motherParticle.GetE();
  m_calculated_mother_p = motherParticle.GetP();
  m_calculated_mother_p_err = motherParticle.GetErrP();
  m_calculated_mother_pt = motherParticle.GetPt();
  m_calculated_mother_pt_err = motherParticle.GetErrPt();
  m_calculated_mother_q = motherParticle.Q();
  // m_calculated_mother_eta = motherParticle.GetEta();
  motherParticle.GetEta(m_calculated_mother_eta, TempError);
  m_calculated_mother_rapidity = motherParticle.GetRapidity();
  m_calculated_mother_theta = motherParticle.GetTheta();
  // m_calculated_mother_phi = motherParticle.GetPhi();
  motherParticle.GetPhi(m_calculated_mother_phi, TempError);
  m_calculated_mother_v = KFParticle_util::calculateEllipsoidVolume(motherParticle);
  m_calculated_mother_pdgID = motherParticle.GetPDG();
  // m_calculated_mother_cov          = &motherParticle.CovarianceMatrix()[0];
  for (int j = 0; j < 21; ++j)
  {
    m_calculated_mother_cov[j] = motherParticle.GetCovariance(j);
  }

  motherParticle.GetMass(m_calculated_mother_mass, m_calculated_mother_mass_err);

  KFParticle* intermediateArray = intermediates.data();
  if (m_config.hasIntermediates)
  {
    for (int i = 0; i < m_config.numIntermediateStates; ++i)
    {
      m_calculated_intermediate_dira[i] = KFParticle_util::eventDIRA(intermediateArray[i], motherParticle);
      m_calculated_intermediate_fdchi2[i] = KFParticle_util::flightDistanceChi2(intermediateArray[i], motherParticle);
      if (m_config.constrainToPrimaryVertex)
      {
        m_calculated_intermediate_PV_dca[i] = intermediateArray[i].GetDistanceFromVertex(vertex_fillbranch);
        m_calculated_intermediate_PV_dca_sig[i] = intermediateArray[i].GetDeviationFromVertex(vertex_fillbranch);
        m_calculated_intermediate_PV_dca_err[i] = m_calculated_intermediate_PV_dca[i] / std::sqrt(m_calculated_intermediate_PV_dca_sig[i]);
        m_calculated_intermediate_PV_dca_xy[i] = intermediateArray[i].GetDistanceFromVertexXY(vertex_fillbranch);
      }
      m_calculated_intermediate_x[i] = intermediateArray[i].GetX();
      m_calculated_intermediate_y[i] = intermediateArray[i].GetY();
      m_calculated_intermediate_z[i] = intermediateArray[i].GetZ();
      m_calculated_intermediate_px[i] = intermediateArray[i].GetPx();
      m_calculated_intermediate_py[i] = intermediateArray[i].GetPy();
      m_calculated_intermediate_pz[i] = intermediateArray[i].GetPz();
      m_calculated_intermediate_pe[i] = intermediateArray[i].GetE();
      m_calculated_intermediate_p[i] = intermediateArray[i].GetP();
      m_calculated_intermediate_p_err[i] = intermediateArray[i].GetErrP();
      m_calculated_intermediate_pt[i] = intermediateArray[i].GetPt();
      m_calculated_intermediate_pt_err[i] = intermediateArray[i].GetErrPt();
      m_calculated_intermediate_q[i] = intermediateArray[i].Q();  // I used to cast this as an int. clang-tidy want Uchar_t to Char_t to int
      // m_calculated_intermediate_eta[i] = intermediateArray[i].GetEta();
      intermediateArray[i].GetEta(m_calculated_intermediate_eta[i], TempError);
      m_calculated_intermediate_rapidity[i] = intermediateArray[i].GetRapidity();
      m_calculated_intermediate_theta[i] = intermediateArray[i].GetTheta();
      // m_calculated_intermediate_phi[i] = intermediateArray[i].GetPhi();
      intermediateArray[i].GetPhi(m_calculated_intermediate_phi[i], TempError);
      m_calculated_intermediate_v[i] = KFParticle_util::calculateEllipsoidVolume(intermediateArray[i]);
      m_calculated_intermediate_pdgID[i] = intermediateArray[i].GetPDG();
      // m_calculated_intermediate_cov[i]          = &intermediateArray[i].CovarianceMatrix()[0];
      for (int j = 0; j < 21; ++j)
      {
        m_calculated_intermediate_cov[i][j] = intermediateArray[i].GetCovariance(j);
      }

      intermediateArray[i].GetMass(m_calculated_intermediate_mass[i], m_calculated_intermediate_mass_err[i]);
      intermediateArray[i].SetProductionVertex(motherParticle);
      intermediateArray[i].GetLifeTime(m_calculated_intermediate_decaytime[i], m_calculated_intermediate_decaytime_err[i]);
      intermediateArray[i].GetDecayLength(m_calculated_intermediate_decaylength[i], m_calculated_intermediate_decaylength_err[i]);
      intermediateArray[i].GetDecayLengthXY(m_calculated_intermediate_decaylength_xy[i], m_calculated_intermediate_decaylength_xy_err[i]);

      m_calculated_intermediate_decaytime[i] /= KFParticle_util::speedOfLight;
      m_calculated_intermediate_decaytime_err[i] /= KFParticle_util::speedOfLight;
    }
  }

  isTrackEMCalmatch = true;
  bool bunchCrossingisZero = true;
  for (int i = 0; i < m_config.numTracks; ++i)
  {
    m_calculated_daughter_mass[i] = daughterArray[i].GetMass();
    if (m_config.constrainToPrimaryVertex)
    {
      m_calculated_daughter_PV_dca[i] = daughterArray[i].GetDistanceFromVertex(vertex_fillbranch);
      m_calculated_daughter_PV_dca_sig[i] = daughterArray[i].GetDeviationFromVertex(vertex_fillbranch);
      m_calculated_daughter_PV_dca_err[i] = m_calculated_daughter_PV_dca[i] / std::sqrt(m_calculated_daughter_PV_dca_sig[i]);
      m_calculated_daughter_PV_dca_xy[i] = daughterArray[i].GetDistanceFromVertexXY(vertex_fillbranch);
      m_calculated_daughter_PV_dca_xy_sig[i] = daughterArray[i].GetDeviationFromVertexXY(vertex_fillbranch);
    }
    m_calculated_daughter_x[i] = daughterArray[i].GetX();
    m_calculated_daughter_y[i] = daughterArray[i].GetY();
    m_calculated_daughter_z[i] = daughterArray[i].GetZ();
    m_calculated_daughter_px[i] = daughterArray[i].GetPx();
    m_calculated_daughter_py[i] = daughterArray[i].GetPy();
    m_calculated_daughter_pz[i] = daughterArray[i].GetPz();
    m_calculated_daughter_pe[i] = daughterArray[i].GetE();
    m_calculated_daughter_p[i] = daughterArray[i].GetP();
    m_calculated_daughter_p_err[i] = daughterArray[i].GetErrP();
    m_calculated_daughter_pt[i] = daughterArray[i].GetPt();
    m_calculated_daughter_pt_err[i] = daughterArray[i].GetErrPt();
    m_calculated_daughter_q[i] = daughterArray[i].Q();
    // m_calculated_daughter_eta[i] = daughterArray[i].GetEta();
    daughterArray[i].GetEta(m_calculated_daughter_eta[i], TempError);
    m_calculated_daughter_rapidity[i] = daughterArray[i].GetRapidity();
    m_calculated_daughter_theta[i] = daughterArray[i].GetTheta();
    // m_calculated_daughter_phi[i] = daughterArray[i].GetPhi();
    daughterArray[i].GetPhi(m_calculated_daughter_phi[i], TempError);
    m_calculated_daughter_chi2[i] = daughterArray[i].GetChi2();
    m_calculated_daughter_ndof[i] = daughterArray[i].GetNDF();
    m_calculated_daughter_trid[i] = daughterArray[i].Id();
    m_calculated_daughter_pdgID[i] = daughterArray[i].GetPDG();
    // m_calculated_daughter_cov[i]          = &daughterArray[i].CovarianceMatrix()[0];
    for (int j = 0; j < 21; ++j)
    {
      m_calculated_daughter_cov[i][j] = daughterArray[i].GetCovariance(j);
    }

    // Now get bunch crossing number for the daughter particle
    SvtxTrackMap* thisTrackMap = findNode::getClass<SvtxTrackMap>(topNode, m_config.trackMapNodeName);
    SvtxTrack* thisTrack = KFParticle_util::getTrack(daughterArray[i].Id(), thisTrackMap);
    m_calculated_daughter_bunch_crossing[i] = thisTrack->get_crossing();
    if (m_calculated_daughter_bunch_crossing[i] != 0)
    {
      bunchCrossingisZero = false;
    }
    if (m_config.dEdxInfo)
    {
      m_calculated_daughter_dedx[i] = KFParticle_util::get_dEdx(topNode, daughterArray[i], m_config.trackMapNodeName);  // m_config.dEdxInfo defaults to false; run get_dEdx_info() to change this
    }
    // m_calculated_daughter_expected_dedx_pion[i] = KFParticle_Tools::get_dEdx_fitValue((Int_t) daughterArray[i].GetQ() * daughterArray[i].GetP(), 211);
    // m_calculated_daughter_expected_dedx_kaon[i] = KFParticle_Tools::get_dEdx_fitValue((Int_t) daughterArray[i].GetQ() * daughterArray[i].GetP(), 321);
    // m_calculated_daughter_expected_dedx_proton[i] = KFParticle_Tools::get_dEdx_fitValue((Int_t) daughterArray[i].GetQ() * daughterArray[i].GetP(), 2212);

    bool tempEMCalmatch = false;
    if (m_config.caloInfo)
    {
      m_calo.fill(topNode, daughterArray[i], i, tempEMCalmatch, vertex_fillbranch);
      if (m_config.requireTrackEMCalMatch && !tempEMCalmatch)
      {
        isTrackEMCalmatch = false;
      }
    }
    if (m_config.truthMatching)
    {
      m_truth.fill(topNode, daughterArray[i], i, vertex_fillbranch, m_config.constrainToPrimaryVertex);
      m_truth.getHepMCInfo(topNode, daughterArray[i], i);
    }
    if (m_config.detectorInfo)
    {
      m_detector.fill(topNode, daughterArray[i], i);
    }
  }

  KFVertex motherDecayVertex;

  int iter = 0;
  // Calcualte jT wrt their own mother, not grandmother
  for (int k = 0; k < m_config.numIntermediateStates; ++k)
  {
    KFVertex intermediateDecayVertex;
    for (int j = 0; j < m_config.numTracksFromIntermediate[k]; ++j)
    {
      m_calculated_daughter_jt[iter] = calculateJT(intermediateArray[k], daughterArray[iter]);
      intermediateDecayVertex += daughterArray[iter];
      ++iter;
    }
    m_calculated_intermediate_chi2[k] = intermediateDecayVertex.GetChi2();
    m_calculated_intermediate_ndof[k] = intermediateDecayVertex.GetNDF();
    motherDecayVertex += intermediateArray[k];
  }
  for (int k = 0; k < num_remaining_tracks; k++)
  {
    m_calculated_daughter_jt[iter] = calculateJT(motherParticle, daughterArray[iter]);
    motherDecayVertex += daughterArray[iter];
    ++iter;
  }

  m_calculated_mother_chi2 = motherDecayVertex.GetChi2();
  m_calculated_mother_ndof = motherDecayVertex.GetNDF();

  iter = 0;
  for (int i = 0; i < m_config.numTracks; ++i)
  {
    for (int j = i + 1; j < m_config.numTracks; ++j)
    {
      m_daughter_dca[iter] = daughterArray[i].GetDistanceFromParticle(daughterArray[j]);
      m_daughter_dca_xy[iter] = daughterArray[i].GetDistanceFromParticleXY(daughterArray[j]);
      m_daughter_dca_sig[iter] = daughterArray[i].GetDeviationFromParticle(daughterArray[j]);
      m_daughter_dca_sig_xy[iter] = daughterArray[i].GetDeviationFromParticleXY(daughterArray[j]);
      ++iter;
    }
  }

  if (m_config.allPVInfo)
  {
    std::vector<KFParticle> sortedDaughterVector;
    sortedDaughterVector.reserve(m_config.numTracks);
    for (int i = 0; i < m_config.numTracks; ++i)
    {
      sortedDaughterVector.push_back(daughterArray[i]);
    }
    allPVInfo(topNode, motherParticle, sortedDaughterVector, intermediates);
  }

  if (m_config.constrainToPrimaryVertex)
  {
    motherParticle.SetProductionVertex(vertex_fillbranch);
    motherParticle.GetLifeTime(m_calculated_mother_decaytime, m_calculated_mother_decaytime_err);
    motherParticle.GetDecayLength(m_calculated_mother_decaylength, m_calculated_mother_decaylength_err);
    motherParticle.GetDecayLengthXY(m_calculated_mother_decaylength_xy, m_calculated_mother_decaylength_xy_err);

    m_calculated_mother_decaytime /= KFParticle_util::speedOfLight;
    m_calculated_mother_decaytime_err /= KFParticle_util::speedOfLight;

    m_calculated_vertex_x = vertex_fillbranch.GetX();
    m_calculated_vertex_y = vertex_fillbranch.GetY();
    m_calculated_vertex_z = vertex_fillbranch.GetZ();
    m_calculated_vertex_v = KFParticle_util::calculateEllipsoidVolume(vertex_fillbranch);
    m_calculated_vertex_chi2 = vertex_fillbranch.GetChi2();
    m_calculated_vertex_ndof = vertex_fillbranch.GetNDF();

    // it only makes sense to calculate PVID for non-fake vertex
    // (this otherwise crashes if m_config.useFakePrimaryVertex is true in an event with no real vertices)
    m_calculated_vertex_ID = m_config.useFakePrimaryVertex ? -100 : KFParticle_util::getPVID(topNode, vertex_fillbranch, m_config.dontUseGlobalVertex,
                                                                                     m_config.useMbdVertex, m_config.vertexMapNodeName);  // -100 is the error value returned by getPVID
    // m_calculated_vertex_cov          = &vertex_fillbranch.CovarianceMatrix()[0];
    for (int j = 0; j < 6; ++j)
    {
      m_calculated_vertex_cov[j] = vertex_fillbranch.GetCovariance(j);
    }
  }

  m_sv_mass = calc_secondary_vertex_mass_noPID(daughters);

  getTracksFromBC(topNode, m_calculated_daughter_bunch_crossing[0], m_config.vertexMapNodeName, m_multiplicity, m_nPVs);
  m_ncharged_siseed_multiplicity = getNchargedSiSeedMultiplicity(topNode, m_calculated_daughter_bunch_crossing[0], m_config.verbosity);
  // cannot retrieve vertex map info from fake PV, hence the second condition
  if (m_config.constrainToPrimaryVertex && !m_config.useFakePrimaryVertex)
  {
    m_nTracksOfVertex = getTracksFromVertex(topNode, vertex_fillbranch, m_config.vertexMapNodeName);
  }
  else
  {
    m_nTracksOfVertex = 0;
  }

  PHNodeIterator nodeIter(topNode);

  PHNode* evtNode = nodeIter.findFirst("EventHeader");

  m_runNumber = -1;
  m_evtNumber = -1;
  m_bco = -1;
  if (evtNode)
  {
    EventHeader* evtHeader = findNode::getClass<EventHeader>(topNode, "EventHeader");
    if (evtHeader)
    {
      m_runNumber = evtHeader->get_RunNumber();
      m_evtNumber = evtHeader->get_EvtSequence();
    }

    auto* gl1packet = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
    if (!gl1packet)
    {
      gl1packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    }

    if (gl1packet)
    {
      m_bco = gl1packet->lValue(0, "BCO") + m_calculated_daughter_bunch_crossing[0];
    }
  }

  if (m_config.useCentrality)
  {
    m_CentInfo = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
    if (!m_CentInfo)
    {
      m_centrality_mbd = -9999.;
    }
    else if (!m_CentInfo->has_centrality_bin(CentralityInfo::PROP::mbd_NS))
    {
      if (m_config.verbosity > 10)
      {
        std::cout << "KFParticle - No centrality information found in CentralityInfo. Setting centrality_MBD to -999." << std::endl;
        m_CentInfo->identify();
      }
      m_centrality_mbd = -999.;
    }
    else if (!bunchCrossingisZero)
    {
      m_centrality_mbd = std::numeric_limits<float>::quiet_NaN();
    }
    else
    {
      m_centrality_mbd = m_CentInfo->get_centrality_bin(CentralityInfo::PROP::mbd_NS);
    }
  }

  if (m_trigger.isAvailable())
  {
    m_trigger.fill(topNode);
  }

  if (fillConditionMet())
  {
    m_tree->Fill();
  }

  // Each helper empties the per-candidate vectors it filled
  if (m_config.caloInfo)
  {
    m_calo.clear(m_config.numTracks);
  }
  if (m_config.truthMatching)
  {
    m_truth.clear(m_config.numTracks);
  }
  if (m_config.detectorInfo)
  {
    m_detector.clear(m_config.numTracks);
  }
  if (m_config.allPVInfo)
  {
    clearAllPVInfo();
  }

  if (m_trigger.isAvailable())
  {
    m_trigger.clear();
  }
}


void KFParticle_nTuple::initializeMomentumCalibrationBranches(const std::string& mother_name)
{
  m_tree->Branch(TString(mother_name) + "_mass", &m_calculated_mother_mass, TString(mother_name) + "_mass/F");
  m_tree->Branch(TString(mother_name) + "_pT", &m_calculated_mother_pt, TString(mother_name) + "_pT/F");
  m_tree->Branch(TString(mother_name) + "_pseudorapidity", &m_calculated_mother_eta, TString(mother_name) + "_pseudorapidity/F");
  m_tree->Branch(TString(mother_name) + "_phi", &m_calculated_mother_phi, TString(mother_name) + "_phi/F");
  m_tree->Branch(TString(mother_name) + "_PDG_ID", &m_calculated_mother_pdgID, TString(mother_name) + "_PDG_ID/I");

  for (int i = 0; i < m_config.numTracks; ++i)
  {
    std::string daughter_number = "track_" + std::to_string(i + 1);

    m_tree->Branch(TString(daughter_number) + "_pT", &m_calculated_daughter_pt[i], TString(daughter_number) + "_pT/F");
    m_tree->Branch(TString(daughter_number) + "_pTErr", &m_calculated_daughter_pt_err[i], TString(daughter_number) + "_pTErr/F");
    m_tree->Branch(TString(daughter_number) + "_pseudorapidity", &m_calculated_daughter_eta[i], TString(daughter_number) + "_pseudorapidity/F");
    m_tree->Branch(TString(daughter_number) + "_phi", &m_calculated_daughter_phi[i], TString(daughter_number) + "_phi/F");
    m_tree->Branch(TString(daughter_number) + "_charge", &m_calculated_daughter_q[i], TString(daughter_number) + "_charge/B");
    m_tree->Branch(TString(daughter_number) + "_PDG_ID", &m_calculated_daughter_pdgID[i], TString(daughter_number) + "_PDG_ID/I");
    m_tree->Branch(TString(daughter_number) + "_Covariance", &m_calculated_daughter_cov[i], TString(daughter_number) + "_Covariance[21]/F", 21);
  }

  m_tree->Branch("runNumber", &m_runNumber, "runNumber/I");
}

void KFParticle_nTuple::fillMomentumCalibrationBranch(PHCompositeNode* topNode,
                                                      KFParticle& motherParticle,
                                                      KFParticle* daughterArray)
{
  motherParticle.GetMass(m_calculated_mother_mass, m_calculated_mother_mass_err);
  m_calculated_mother_pt = motherParticle.GetPt();
  motherParticle.GetEta(m_calculated_mother_eta, TempError);
  motherParticle.GetPhi(m_calculated_mother_phi, TempError);
  m_calculated_mother_pdgID = motherParticle.GetPDG();

  for (int i = 0; i < m_config.numTracks; ++i)
  {
    m_calculated_daughter_pt[i] = daughterArray[i].GetPt();
    m_calculated_daughter_pt_err[i] = daughterArray[i].GetErrPt();
    daughterArray[i].GetEta(m_calculated_daughter_eta[i], TempError);
    daughterArray[i].GetPhi(m_calculated_daughter_phi[i], TempError);
    m_calculated_daughter_q[i] = daughterArray[i].Q();
    m_calculated_daughter_pdgID[i] = daughterArray[i].GetPDG();
    for (int j = 0; j < 21; ++j)
    {
      m_calculated_daughter_cov[i][j] = daughterArray[i].GetCovariance(j);
    }
  }

  EventHeader* evtHeader = findNode::getClass<EventHeader>(topNode, "EventHeader");
  m_runNumber = evtHeader ? evtHeader->get_RunNumber() : -1;

  m_tree->Fill();
}

float KFParticle_nTuple::calc_secondary_vertex_mass_noPID(std::vector<KFParticle> kfp_daughters) const
{
  KFParticle mother_noPID;
  KFParticle* daughterArray = kfp_daughters.data();

  for (int i = 0; i < m_config.numTracks; ++i)
  {
    KFParticle daughter_noPID;
    float f_trackParameters[6];
    float f_trackCovariance[21];
    for (int j = 0; j < 6; ++j)
    {
      f_trackParameters[j] = daughterArray[i].GetParameter(j);
    }
    for (int j = 0; j < 21; ++j)
    {
      f_trackCovariance[j] = daughterArray[i].GetCovariance(j);
    }
    daughter_noPID.Create(f_trackParameters, f_trackCovariance, daughterArray[i].Q(), -1);
    mother_noPID.AddDaughter(daughter_noPID);
  }

  return mother_noPID.GetMass();
}

bool KFParticle_nTuple::fillConditionMet() const
{
  // Fill if this is a track-only analysis, if track-calo matching is not required,
  // or otherwise only if the track was matched to an EMCal cluster
  return !m_config.caloInfo || !m_config.requireTrackEMCalMatch || isTrackEMCalmatch;
}

void KFParticle_nTuple::allPVInfo(PHCompositeNode *topNode,
                                  const KFParticle &motherParticle,
                                  std::vector<KFParticle> daughters,
                                  std::vector<KFParticle> intermediates)
{
  std::vector<KFParticle> primaryVertices = KFParticle_util::makeAllPrimaryVertices(topNode, m_config);

  for (auto &primaryVertice : primaryVertices)
  {
    allPV_x.push_back(primaryVertice.GetX());
    allPV_y.push_back(primaryVertice.GetY());
    allPV_z.push_back(primaryVertice.GetZ());

    allPV_mother_PV_DCA.push_back(motherParticle.GetDistanceFromVertex(primaryVertice));
    allPV_mother_PV_DCA_StdDev.push_back(motherParticle.GetDeviationFromVertex(primaryVertice));

    for (unsigned int j = 0; j < daughters.size(); ++j)
    {
      allPV_daughter_PV_DCA[j].push_back(daughters[j].GetDistanceFromVertex(primaryVertice));
      allPV_daughter_PV_DCA_StdDev[j].push_back(daughters[j].GetDeviationFromVertex(primaryVertice));
    }

    for (unsigned int j = 0; j < intermediates.size(); ++j)
    {
      allPV_intermediates_PV_DCA[j].push_back(intermediates[j].GetDistanceFromVertex(primaryVertice));
      allPV_intermediates_PV_DCA_StdDev[j].push_back(intermediates[j].GetDeviationFromVertex(primaryVertice));
    }
  }
}

void KFParticle_nTuple::clearAllPVInfo()
{
  for (int i = 0; i < m_config.numTracks; ++i)
  {
    allPV_daughter_PV_DCA[i].clear();
    allPV_daughter_PV_DCA_StdDev[i].clear();
  }

  allPV_x.clear();
  allPV_y.clear();
  allPV_z.clear();

  allPV_mother_PV_DCA.clear();
  allPV_mother_PV_DCA_StdDev.clear();

  for (int i = 0; i < m_config.numIntermediateStates; ++i)
  {
    allPV_intermediates_PV_DCA[i].clear();
    allPV_intermediates_PV_DCA_StdDev[i].clear();
  }
}

void KFParticle_nTuple::getTracksFromBC(PHCompositeNode *topNode, const int &bunch_crossing, const std::string &vertexMapName, int &nTracks, int &nPVs) const
{
  if (m_config.useMbdVertex)  // If you're using the MBD vertex then there is no way to know which tracks are associated to it
  {
    return;
  }

  const std::string &vtxMN = vertexMapName.empty() ? m_config.vertexMapNodeName : vertexMapName;

  nTracks = 0;
  nPVs = 0;
  auto *dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, vtxMN);
  for (SvtxVertexMap::ConstIter iter = dst_vertexmap->begin(); iter != dst_vertexmap->end(); ++iter)
  {
    SvtxVertex *dst_vertex = iter->second;
    if ((int) dst_vertex->get_beam_crossing() == bunch_crossing)
    {
      nTracks += dst_vertex->size_tracks();
      ++nPVs;
    }
  }
}

int KFParticle_nTuple::getTracksFromVertex(PHCompositeNode *topNode, const KFParticle &vertex, const std::string &vertexMapName) const
{
  if (m_config.useMbdVertex)  // If you're using the MBD vertex then there is no way to know which tracks are associated to it
  {
    return 0;
  }

  const std::string &vtxMN = vertexMapName.empty() ? m_config.vertexMapNodeName : vertexMapName;

  SvtxVertex *associated_vertex = nullptr;
  if (m_config.dontUseGlobalVertex)
  {
    auto *dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, vtxMN);
    associated_vertex = dst_vertexmap ? dst_vertexmap->get(vertex.Id()) : nullptr;
  }
  else
  {
    auto *dst_globalvertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    auto *dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, vtxMN);
    if (!dst_globalvertexmap || !dst_vertexmap)
    {
      std::cout << "Can't continue in KFParticle_nTuple::getTracksFromVertex" << std::endl;
      return 0;
    }
    auto *associated_gvertex = dst_globalvertexmap->get(vertex.Id());
    if (!associated_gvertex)
    {
      return 0;
    }

    auto svtxiter = associated_gvertex->find_vertexes(GlobalVertex::SVTX);
    if (svtxiter == associated_gvertex->end_vertexes())
    {
      return 0;
    }
    const auto &svtxvertexvector = svtxiter->second;

    for (const auto &gvertex : svtxvertexvector)
    {
      auto vtxiter = dst_vertexmap->find(gvertex->get_id());
      if (vtxiter != dst_vertexmap->end())
      {
        associated_vertex = vtxiter->second;
      }
    }
  }
  if (associated_vertex)
  {
    return associated_vertex->size_tracks();
  }

  return 0;
}

int KFParticle_nTuple::getNchargedSiSeedMultiplicity(PHCompositeNode *topNode, const int &bunch_crossing, int verbosity)
{
  auto *m_siliconSeeds = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  //auto clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_SEED");
  if (!m_siliconSeeds)
  {
    if (verbosity > 0)
    {
      std::cout << PHWHERE << " ERROR: Can't find SiliconTrackSeedContainer " << std::endl;
    }
    return -1;
  }
  //if (!clustermap)
  //{
  //  std::cout << " ERROR: Can't find TRKR_CLUSTER " << std::endl;
  //  return -1;
  //}

  TrackSeed *_tracklet_si;
  int ncharged_multiplicity = 0;
  const int nMapsCut = 1;
  const int nInttCut = 1;
  for (unsigned int phtrk_iter_si = 0;
         phtrk_iter_si < m_siliconSeeds->size();
         ++phtrk_iter_si)
  {
    _tracklet_si = m_siliconSeeds->get(phtrk_iter_si);
    if (!_tracklet_si)
    {
      continue;
    }
    if (_tracklet_si->get_crossing() == bunch_crossing)
    {
      if (TrackAnalysisUtils::get_cluster_count(_tracklet_si,TrkrDefs::mvtxId) >= nMapsCut &&
          TrackAnalysisUtils::get_cluster_count(_tracklet_si,TrkrDefs::inttId) >= nInttCut)
      {
        ncharged_multiplicity++;
      }
    }
  }

  return ncharged_multiplicity;
}

float KFParticle_nTuple::calculateJT(const KFParticle &mother, const KFParticle &daughter)
{
  Eigen::Vector3f motherP = Eigen::Vector3f(mother.GetPx(), mother.GetPy(), mother.GetPz());
  Eigen::Vector3f daughterP = Eigen::Vector3f(daughter.GetPx(), daughter.GetPy(), daughter.GetPz());

  Eigen::Vector3f motherP_X_daughterP = motherP.cross(daughterP);

  float jT = (motherP_X_daughterP.norm()) / motherP.norm();

  return jT;
}
