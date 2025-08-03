#include "KFParticle_nTuple.h"

#include "KFParticle_Tools.h"

#include <ffaobjects/EventHeader.h>
#include <ffarawobjects/Gl1Packet.h>

#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>
#include <trackbase_historic/SvtxTrack.h>   
#include <trackbase_historic/SvtxTrackMap.h>

#include <KFParticle.h>
#include <KFVertex.h>

#include <Rtypes.h>
#include <TString.h>  // for TString, operator+
#include <TTree.h>

#include <algorithm>  // for max
#include <cmath>
#include <cstdlib>  // for abs, size_t
#include <map>      // for map, _Rb_tree_iterator, map<>:...

class PHCompositeNode;
class PHNode;

/// Create necessary objects
KFParticle_Tools kfpTupleTools;
float TempError;

void KFParticle_nTuple::initializeVariables()
{
  // m_calculated_mother_cov = -99;
}

void KFParticle_nTuple::initializeBranches(PHCompositeNode* topNode)
{
  //kfpTupleTools.init_dEdx_fits(); //Cant do this! Two trees open at once!

  delete m_tree;
  m_tree = new TTree("DecayTree", "DecayTree");
  m_tree->OptimizeBaskets();
  m_tree->SetAutoSave(-5e6);  // Save the output file every 5MB

  std::string mother_name;
  if (m_mother_name.empty())
  {
    mother_name = "mother";
  }
  else
  {
    mother_name = m_mother_name;
  }

  size_t pos;
  std::string undrscr = "_";
  std::string nothing = "";
  std::map<std::string, std::string> forbiddenStrings;
  forbiddenStrings["/"] = undrscr;
  forbiddenStrings["("] = undrscr;
  forbiddenStrings[")"] = nothing;
  forbiddenStrings["+"] = "plus";
  forbiddenStrings["-"] = "minus";
  forbiddenStrings["*"] = "star";
  for (auto const& [badString, goodString] : forbiddenStrings)
  {
    while ((pos = mother_name.find(badString)) != std::string::npos)
    {
      mother_name.replace(pos, 1, goodString);
    }
  }

  m_tree->Branch(TString(mother_name) + "_mass", &m_calculated_mother_mass, TString(mother_name) + "_mass/F");
  m_tree->Branch(TString(mother_name) + "_massErr", &m_calculated_mother_mass_err, TString(mother_name) + "_massErr/F");
  if (m_constrain_to_vertex_nTuple)
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
    m_tree->Branch(TString(mother_name) + "_IP", &m_calculated_mother_ip, TString(mother_name) + "_IP/F");
    m_tree->Branch(TString(mother_name) + "_IPchi2", &m_calculated_mother_ipchi2, TString(mother_name) + "_IPchi2/F");
    m_tree->Branch(TString(mother_name) + "_IPErr", &m_calculated_mother_ip_err, TString(mother_name) + "_IPErr/F");
    m_tree->Branch(TString(mother_name) + "_IP_xy", &m_calculated_mother_ip_xy, TString(mother_name) + "_IP_xy/F");
  }
  if (m_get_all_PVs)
  {
    m_tree->Branch(TString(mother_name) + "_IP_allPV", &allPV_mother_IP);
    m_tree->Branch(TString(mother_name) + "_IPchi2_allPV", &allPV_mother_IPchi2);
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
  if (m_has_intermediates_nTuple)
  {
    for (int i = 0; i < m_num_intermediate_states_nTuple; ++i)
    {
      std::string intermediate_name = m_intermediate_name_ntuple[i];

      // Note, TBranch will not allow the leaf to contain a forward slash as it is used to define the branch type. Causes problems with J/psi
      for (auto const& [badString, goodString] : forbiddenStrings)
      {
        while ((pos = intermediate_name.find(badString)) != std::string::npos)
        {
          intermediate_name.replace(pos, 1, goodString);
        }
      }

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
      if (m_constrain_to_vertex_nTuple)
      {
        m_tree->Branch(TString(intermediate_name) + "_IP", &m_calculated_intermediate_ip[i], TString(intermediate_name) + "_IP/F");
        m_tree->Branch(TString(intermediate_name) + "_IPchi2", &m_calculated_intermediate_ipchi2[i], TString(intermediate_name) + "_IPchi2/F");
        m_tree->Branch(TString(intermediate_name) + "_IPErr", &m_calculated_intermediate_ip_err[i], TString(intermediate_name) + "_IPErr/F");
        m_tree->Branch(TString(intermediate_name) + "_IP_xy", &m_calculated_intermediate_ip_xy[i], TString(intermediate_name) + "_IP_xy/F");
      }
      if (m_get_all_PVs)
      {
        m_tree->Branch(TString(intermediate_name) + "_IP_allPV", &allPV_intermediates_IP[i]);
        m_tree->Branch(TString(intermediate_name) + "_IPchi2_allPV", &allPV_intermediates_IPchi2[i]);
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

      for (int j = 0; j < m_num_tracks_from_intermediate_nTuple[i]; ++j)
      {
        intermediateNameMapping.push_back(intermediate_name + "_");
      }
    }
  }

  int num_intermediate_tracks = 0;
  for (int i = 0; i < m_num_intermediate_states_nTuple; ++i)
  {
    num_intermediate_tracks += m_num_tracks_from_intermediate_nTuple[i];
  }

  for (int i = 0; i < m_num_tracks_nTuple; ++i)
  {
    std::string daughter_number = "track_" + std::to_string(i + 1);

    if (m_has_intermediates_nTuple && i < num_intermediate_tracks)
    {
      daughter_number.insert(0, intermediateNameMapping[i]);
    }

    m_tree->Branch(TString(daughter_number) + "_mass", &m_calculated_daughter_mass[i], TString(daughter_number) + "_mass/F");
    if (m_constrain_to_vertex_nTuple)
    {
      m_tree->Branch(TString(daughter_number) + "_IP", &m_calculated_daughter_ip[i], TString(daughter_number) + "_IP/F");
      m_tree->Branch(TString(daughter_number) + "_IPchi2", &m_calculated_daughter_ipchi2[i], TString(daughter_number) + "_IPchi2/F");
      m_tree->Branch(TString(daughter_number) + "_IPErr", &m_calculated_daughter_ip_err[i], TString(daughter_number) + "_IPErr/F");
      m_tree->Branch(TString(daughter_number) + "_IP_xy", &m_calculated_daughter_ip_xy[i], TString(daughter_number) + "_IP_xy/F");
    }
    if (m_get_all_PVs)
    {
      m_tree->Branch(TString(daughter_number) + "_IP_allPV", &allPV_daughter_IP[i]);
      m_tree->Branch(TString(daughter_number) + "_IPchi2_allPV", &allPV_daughter_IPchi2[i]);
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
    m_tree->Branch(TString(daughter_number) + "_calculated_dEdx", &m_calculated_daughter_dedx[i], TString(daughter_number) + "_calculated_dEdx/F");
    //m_tree->Branch(TString(daughter_number) + "_expected_pion_dEdx", &m_calculated_daughter_expected_dedx_pion[i], TString(daughter_number) + "_expected_pion_dEdx/F");
    //m_tree->Branch(TString(daughter_number) + "_expected_kaon_dEdx", &m_calculated_daughter_expected_dedx_kaon[i], TString(daughter_number) + "_expected_kaon_dEdx/F");
    //m_tree->Branch(TString(daughter_number) + "_expected_proton_dEdx", &m_calculated_daughter_expected_dedx_proton[i], TString(daughter_number) + "_expected_proton_dEdx/F");

    if (m_calo_info)
    {
      initializeCaloBranches(m_tree, i, daughter_number);
    }
    if (m_truth_matching)
    {
      initializeTruthBranches(m_tree, i, daughter_number, m_constrain_to_vertex_nTuple);
    }
    if (m_detector_info)
    {
      initializeDetectorBranches(m_tree, i, daughter_number);
    }
  }

  int iter = 0;
  for (int i = 0; i < m_num_tracks_nTuple; ++i)
  {
    for (int j = 0; j < m_num_tracks_nTuple; ++j)
    {
      if (i < j)
      {
        std::string dca_branch_name = "track_" + std::to_string(i + 1) + "_track_" + std::to_string(j + 1) + "_DCA";
        std::string dca_leaf_name = dca_branch_name + "/F";
        m_tree->Branch(dca_branch_name.c_str(), &m_daughter_dca[iter], dca_leaf_name.c_str());

        std::string dca_branch_name_xy = dca_branch_name + "_xy";
        std::string dca_leaf_name_xy = dca_branch_name_xy + "/F";
        m_tree->Branch(dca_branch_name_xy.c_str(), &m_daughter_dca_xy[iter], dca_leaf_name_xy.c_str());

        ++iter;
      }
    }
  }

  if (m_constrain_to_vertex_nTuple)
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
  if (m_get_all_PVs)
  {
    m_tree->Branch("all_primary_vertex_x", &allPV_x);
    m_tree->Branch("all_primary_vertex_y", &allPV_y);
    m_tree->Branch("all_primary_vertex_z", &allPV_z);
  }

  m_tree->Branch("secondary_vertex_mass_pionPID", &m_sv_mass, "secondary_vertex_mass_pionPID/F");

  m_tree->Branch("nPrimaryVerticesOfBC", &m_nPVs, "nPrimaryVerticesOfBC/I");
  m_tree->Branch("nTracksOfBC", &m_multiplicity, "nTracksOfBC/I");
  m_tree->Branch("nTracksOfVertex", &m_nTracksOfVertex, "nTracksOfVertex/I");

  m_tree->Branch("runNumber", &m_runNumber, "runNumber/I");
  m_tree->Branch("eventNumber", &m_evtNumber, "eventNumber/I");
  m_tree->Branch("BCO", &m_bco, "BCO/L");

  if (m_get_trigger_info)
  {
    m_trigger_info_available = buildTriggerBranches(topNode, m_tree);
  }
}

void KFParticle_nTuple::fillBranch(PHCompositeNode* topNode,
                                   KFParticle motherParticle,
                                   const KFParticle& vertex_fillbranch,
                                   std::vector<KFParticle> daughters,
                                   std::vector<KFParticle> intermediates)
{
  const float speedOfLight = 2.99792458e-2;

  KFParticle temp;
  KFParticle* daughterArray = &daughters[0];

  bool switchTrackPosition;

  int num_tracks_used_by_intermediates = 0;
  for (int k = 0; k < m_num_intermediate_states_nTuple; ++k)  // Rearrange tracks from intermediate states
  {
    for (int i = 0; i < m_num_tracks_from_intermediate_nTuple[k]; ++i)
    {
      for (int j = i + 1; j < m_num_tracks_from_intermediate_nTuple[k]; ++j)
      {
        int particleAElement = i + num_tracks_used_by_intermediates;
        int particleBElement = j + num_tracks_used_by_intermediates;
        int particleAPID = daughterArray[particleAElement].GetPDG();
        int particleBPID = daughterArray[particleBElement].GetPDG();

        if (m_get_charge_conjugate_nTuple)
        {
          float daughterA_mass = kfpTupleTools.getParticleMass(particleAPID);
          float daughterB_mass = kfpTupleTools.getParticleMass(particleBPID);
          switchTrackPosition = daughterA_mass > daughterB_mass;
        }
        else
        {
          switchTrackPosition = particleAPID > particleBPID;
        }
        if (switchTrackPosition)
        {
          temp = daughterArray[particleAElement];
          daughterArray[particleAElement] = daughterArray[particleBElement];
          daughterArray[particleBElement] = temp;
        }
      }
    }
    num_tracks_used_by_intermediates += m_num_tracks_from_intermediate_nTuple[k];
  }

  int num_remaining_tracks = m_num_tracks_nTuple - num_tracks_used_by_intermediates;

  for (int i = 0; i < num_remaining_tracks; i++)
  {
    for (int j = i + 1; j < num_remaining_tracks; j++)
    {
      int particleAElement = i + num_tracks_used_by_intermediates;
      int particleBElement = j + num_tracks_used_by_intermediates;
      int particleAPID = daughterArray[particleAElement].GetPDG();
      int particleBPID = daughterArray[particleBElement].GetPDG();

      if (m_get_charge_conjugate_nTuple)
      {
        float daughterA_mass = kfpTupleTools.getParticleMass(particleAPID);
        float daughterB_mass = kfpTupleTools.getParticleMass(particleBPID);
        switchTrackPosition = daughterA_mass > daughterB_mass;
      }
      else
      {
        switchTrackPosition = particleAPID > particleBPID;
      }
      if (switchTrackPosition)
      {
        temp = daughterArray[particleAElement];
        daughterArray[particleAElement] = daughterArray[particleBElement];
        daughterArray[particleBElement] = temp;
      }
    }
  }

  if (m_extrapolateTracksToSV_nTuple)
  {
    for (unsigned int i = 0; i < daughters.size(); ++i)
    {
      daughterArray[i].SetProductionVertex(motherParticle);
    }
  }

  if (m_constrain_to_vertex_nTuple)
  {
    m_calculated_mother_dira = kfpTupleTools.eventDIRA(motherParticle, vertex_fillbranch);
    m_calculated_mother_dira_xy = kfpTupleTools.eventDIRA(motherParticle, vertex_fillbranch, false);
    m_calculated_mother_fdchi2 = kfpTupleTools.flightDistanceChi2(motherParticle, vertex_fillbranch);
    m_calculated_mother_ip = motherParticle.GetDistanceFromVertex(vertex_fillbranch);
    m_calculated_mother_ipchi2 = motherParticle.GetDeviationFromVertex(vertex_fillbranch);
    m_calculated_mother_ip_err = m_calculated_mother_ip / std::sqrt(m_calculated_mother_ipchi2);
    m_calculated_mother_ip_xy = motherParticle.GetDistanceFromVertexXY(vertex_fillbranch);
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
  //m_calculated_mother_eta = motherParticle.GetEta();
  motherParticle.GetEta(m_calculated_mother_eta,TempError);
  m_calculated_mother_rapidity = motherParticle.GetRapidity();
  m_calculated_mother_theta = motherParticle.GetTheta();
  //m_calculated_mother_phi = motherParticle.GetPhi();
  motherParticle.GetPhi(m_calculated_mother_phi, TempError);
  m_calculated_mother_v = kfpTupleTools.calculateEllipsoidVolume(motherParticle);
  m_calculated_mother_pdgID = motherParticle.GetPDG();
  // m_calculated_mother_cov          = &motherParticle.CovarianceMatrix()[0];
  for (int j = 0; j < 21; ++j)
  {
    m_calculated_mother_cov[j] = motherParticle.GetCovariance(j);
  }

  motherParticle.GetMass(m_calculated_mother_mass, m_calculated_mother_mass_err);

  KFParticle* intermediateArray = &intermediates[0];
  if (m_has_intermediates_nTuple)
  {
    for (int i = 0; i < m_num_intermediate_states_nTuple; ++i)
    {
      m_calculated_intermediate_dira[i] = kfpTupleTools.eventDIRA(intermediateArray[i], motherParticle);
      m_calculated_intermediate_fdchi2[i] = kfpTupleTools.flightDistanceChi2(intermediateArray[i], motherParticle);
      if (m_constrain_to_vertex_nTuple)
      {
        m_calculated_intermediate_ip[i] = intermediateArray[i].GetDistanceFromVertex(vertex_fillbranch);
        m_calculated_intermediate_ipchi2[i] = intermediateArray[i].GetDeviationFromVertex(vertex_fillbranch);
        m_calculated_intermediate_ip_err[i] = m_calculated_intermediate_ip[i] / std::sqrt(m_calculated_intermediate_ipchi2[i]);
        m_calculated_intermediate_ip_xy[i] = intermediateArray[i].GetDistanceFromVertexXY(vertex_fillbranch);
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
      //m_calculated_intermediate_eta[i] = intermediateArray[i].GetEta();
      intermediateArray[i].GetEta(m_calculated_intermediate_eta[i], TempError);
      m_calculated_intermediate_rapidity[i] = intermediateArray[i].GetRapidity();
      m_calculated_intermediate_theta[i] = intermediateArray[i].GetTheta();
      //m_calculated_intermediate_phi[i] = intermediateArray[i].GetPhi();
      intermediateArray[i].GetPhi(m_calculated_intermediate_phi[i], TempError);
      m_calculated_intermediate_v[i] = kfpTupleTools.calculateEllipsoidVolume(intermediateArray[i]);
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

      m_calculated_intermediate_decaytime[i] /= speedOfLight;
      m_calculated_intermediate_decaytime_err[i] /= speedOfLight;
    }
  }

  isTrackEMCalmatch = true;
  for (int i = 0; i < m_num_tracks_nTuple; ++i)
  {
    m_calculated_daughter_mass[i] = daughterArray[i].GetMass();
    if (m_constrain_to_vertex_nTuple)
    {
      m_calculated_daughter_ip[i] = daughterArray[i].GetDistanceFromVertex(vertex_fillbranch);
      m_calculated_daughter_ipchi2[i] = daughterArray[i].GetDeviationFromVertex(vertex_fillbranch);
      m_calculated_daughter_ip_err[i] = m_calculated_daughter_ip[i] / std::sqrt(m_calculated_daughter_ipchi2[i]);
      m_calculated_daughter_ip_xy[i] = daughterArray[i].GetDistanceFromVertexXY(vertex_fillbranch);
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
    //m_calculated_daughter_eta[i] = daughterArray[i].GetEta();
    daughterArray[i].GetEta(m_calculated_daughter_eta[i], TempError);
    m_calculated_daughter_rapidity[i] = daughterArray[i].GetRapidity();
    m_calculated_daughter_theta[i] = daughterArray[i].GetTheta();
    //m_calculated_daughter_phi[i] = daughterArray[i].GetPhi();
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

    //Now get bunch crossing number for the daughter particle
    SvtxTrackMap *thisTrackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name_nTuple.c_str());
    SvtxTrack *thisTrack = getTrack(daughterArray[i].Id(), thisTrackMap);
    m_calculated_daughter_bunch_crossing[i] = thisTrack->get_crossing(); 
    m_calculated_daughter_dedx[i] = kfpTupleTools.get_dEdx(topNode, daughterArray[i]);
    //m_calculated_daughter_expected_dedx_pion[i] = kfpTupleTools.get_dEdx_fitValue((Int_t) daughterArray[i].GetQ() * daughterArray[i].GetP(), 211);
    //m_calculated_daughter_expected_dedx_kaon[i] = kfpTupleTools.get_dEdx_fitValue((Int_t) daughterArray[i].GetQ() * daughterArray[i].GetP(), 321);
    //m_calculated_daughter_expected_dedx_proton[i] = kfpTupleTools.get_dEdx_fitValue((Int_t) daughterArray[i].GetQ() * daughterArray[i].GetP(), 2212);

    bool tempEMCalmatch = false;
    if (m_calo_info)
    {
      fillCaloBranch(topNode, m_tree, daughterArray[i], i, tempEMCalmatch);
      if (m_require_track_emcal_match && !tempEMCalmatch)
      {
        isTrackEMCalmatch = false;
      }
    }
    if (m_truth_matching)
    {
      fillTruthBranch(topNode, m_tree, daughterArray[i], i, vertex_fillbranch, m_constrain_to_vertex_nTuple);
    }
    if (m_truth_matching)
    {
      getHepMCInfo(topNode, m_tree, daughterArray[i], i);
    }
    if (m_detector_info)
    {
      fillDetectorBranch(topNode, m_tree, daughterArray[i], i);
    }
  }

  KFVertex motherDecayVertex;

  int iter = 0;
  // Calcualte jT wrt their own mother, not grandmother
  for (int k = 0; k < m_num_intermediate_states_nTuple; ++k)
  {
    KFVertex intermediateDecayVertex;
    for (int j = 0; j < m_num_tracks_from_intermediate_nTuple[k]; ++j)
    {
      m_calculated_daughter_jt[iter] = kfpTupleTools.calculateJT(intermediateArray[k], daughterArray[iter]);
      intermediateDecayVertex += daughterArray[iter];
      ++iter;
    }
    m_calculated_intermediate_chi2[k] = intermediateDecayVertex.GetChi2();
    m_calculated_intermediate_ndof[k] = intermediateDecayVertex.GetNDF();
    motherDecayVertex += intermediateArray[k];
  }
  for (int k = 0; k < num_remaining_tracks; k++)
  {
    m_calculated_daughter_jt[iter] = kfpTupleTools.calculateJT(motherParticle, daughterArray[iter]);
    motherDecayVertex += daughterArray[iter];
    ++iter;
  }

  m_calculated_mother_chi2 = motherDecayVertex.GetChi2();
  m_calculated_mother_ndof = motherDecayVertex.GetNDF();

  iter = 0;
  for (int i = 0; i < m_num_tracks_nTuple; ++i)
  {
    for (int j = 0; j < m_num_tracks_nTuple; ++j)
    {
      if (i < j)
      {
        m_daughter_dca[iter] = daughterArray[i].GetDistanceFromParticle(daughterArray[j]);
        m_daughter_dca_xy[iter] = daughterArray[i].GetDistanceFromParticleXY(daughterArray[j]);
        ++iter;
      }
    }
  }

  if (m_get_all_PVs)
  {
    std::vector<KFParticle> sortedDaughterVector;
    sortedDaughterVector.reserve(m_num_tracks_nTuple);
    for (int i = 0; i < m_num_tracks_nTuple; ++i)
    {
      sortedDaughterVector.push_back(daughterArray[i]);
    }
    allPVInfo(topNode, m_tree, motherParticle, sortedDaughterVector, intermediates);
  }

  if (m_constrain_to_vertex_nTuple)
  {
    motherParticle.SetProductionVertex(vertex_fillbranch);
    motherParticle.GetLifeTime(m_calculated_mother_decaytime, m_calculated_mother_decaytime_err);
    motherParticle.GetDecayLength(m_calculated_mother_decaylength, m_calculated_mother_decaylength_err);
    motherParticle.GetDecayLengthXY(m_calculated_mother_decaylength_xy, m_calculated_mother_decaylength_xy_err);

    m_calculated_mother_decaytime /= speedOfLight;
    m_calculated_mother_decaytime_err /= speedOfLight;

    m_calculated_vertex_x = vertex_fillbranch.GetX();
    m_calculated_vertex_y = vertex_fillbranch.GetY();
    m_calculated_vertex_z = vertex_fillbranch.GetZ();
    m_calculated_vertex_v = kfpTupleTools.calculateEllipsoidVolume(vertex_fillbranch);
    m_calculated_vertex_chi2 = vertex_fillbranch.GetChi2();
    m_calculated_vertex_ndof = vertex_fillbranch.GetNDF();

    // it only makes sense to calculate PVID for non-fake vertex
    // (this otherwise crashes if m_use_fake_pv_nTuple is true in an event with no real vertices)
    if(m_use_fake_pv_nTuple) m_calculated_vertex_ID = -100; // error value returned by getPVID
    else m_calculated_vertex_ID = getPVID(topNode, vertex_fillbranch);
    // m_calculated_vertex_cov          = &vertex_fillbranch.CovarianceMatrix()[0];
    for (int j = 0; j < 6; ++j)
    {
      m_calculated_vertex_cov[j] = vertex_fillbranch.GetCovariance(j);
    }
  }

  m_sv_mass = calc_secondary_vertex_mass_noPID(daughters);

  kfpTupleTools.getTracksFromBC(topNode, m_calculated_daughter_bunch_crossing[0], m_vtx_map_node_name_nTuple, m_multiplicity, m_nPVs);  
  // cannot retrieve vertex map info from fake PV, hence the second condition
  if (m_constrain_to_vertex_nTuple && !m_use_fake_pv_nTuple)
  {
    m_nTracksOfVertex = kfpTupleTools.getTracksFromVertex(topNode, vertex_fillbranch, m_vtx_map_node_name_nTuple);
  }
  else
  {
    m_nTracksOfVertex = 0;
  }

  PHNodeIterator nodeIter(topNode);

  PHNode* evtNode = dynamic_cast<PHNode*>(nodeIter.findFirst("EventHeader"));

  if (evtNode)
  {
    EventHeader* evtHeader = findNode::getClass<EventHeader>(topNode, "EventHeader");
    m_runNumber = evtHeader->get_RunNumber();
    m_evtNumber = evtHeader->get_EvtSequence();

    auto gl1packet = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
    if (!gl1packet)
    {
      gl1packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    }
    m_bco = m_trigger_info_available ? gl1packet->lValue(0, "BCO") + m_calculated_daughter_bunch_crossing[0] : 0;
  }
  else
  {
    m_runNumber = m_evtNumber = m_bco = -1;
  }

  if (m_trigger_info_available)
  {
    fillTriggerBranches(topNode);
  }

  if (fillConditionMet())
  {
    m_tree->Fill();
  }

  if (m_truth_matching || m_detector_info)
  {
    clearVectors();
  }

  if (m_trigger_info_available)
  {
    resetTriggerBranches();
  }
}

float KFParticle_nTuple::calc_secondary_vertex_mass_noPID(std::vector<KFParticle> kfp_daughters)
{
  KFParticle mother_noPID;
  KFParticle* daughterArray = &kfp_daughters[0];

  for (int i = 0; i < m_num_tracks_nTuple; ++i)
  {
    KFParticle daughter_noPID;
    float f_trackParameters[6], f_trackCovariance[21];
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

bool KFParticle_nTuple::fillConditionMet()
{
  // return true if this is a track-only analysis
  if (!m_calo_info)
  {
    return true;
  }

  // return true if do not require track-calo matching
  if (!m_require_track_emcal_match)
  {
    return true;
  }

  // if requiring track-calo matching, the match result is returned
  return isTrackEMCalmatch;
}
