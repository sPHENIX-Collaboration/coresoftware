#ifndef KFPARTICLESPHENIX_KFPARTICLENTUPLE_H
#define KFPARTICLESPHENIX_KFPARTICLENTUPLE_H

#include "KFParticle_truthAndDetTools.h"

#include <KFParticle.h>

#include <string>  // for string
#include <vector>

class PHCompositeNode;
class TTree;

class KFParticle_nTuple : public KFParticle_truthAndDetTools
{
 public:
  ///Constructor
  KFParticle_nTuple();

  ///Destructor
  virtual ~KFParticle_nTuple() {}

  ///Unused for now, variables are initialised in the header
  void initializeVariables();

  ///Initialises required branches based off the user selection (number of tracks, PV constraints etc ) and sets branch names if specified
  void initializeBranches();

  ///Fills required information for your selection, also requests truth and detector information if needed
  void fillBranch(PHCompositeNode *topNode,
                  KFParticle motherParticle,
                  KFParticle vertex,
                  std::vector<KFParticle> daughters,
                  std::vector<KFParticle> intermediates,
                  int nPVs, int multiplicity);

  float calc_secondary_vertex_mass_noPID(std::vector<KFParticle> kfp_daughters);

 protected:
  bool m_has_intermediates_nTuple;
  bool m_constrain_to_vertex_nTuple;
  bool m_get_all_PVs;
  //int m_num_intermediate_states_nTuple;
  //int m_num_tracks_from_intermediate_nTuple[99];
  std::vector<int> m_num_tracks_from_intermediate_nTuple;
  bool m_truth_matching;
  bool m_detector_info;
  bool m_calo_info;
  std::string m_mother_name;
  //std::string m_vtx_map_node_name_nTuple;
  bool m_use_intermediate_name;
  bool m_get_charge_conjugate_nTuple;
  std::vector<std::string> m_intermediate_name_ntuple;

 private:
  TTree *m_tree;

  float m_calculated_mother_mass = -1;
  float m_calculated_mother_mass_err = -1;
  float m_calculated_mother_decaytime = -1;
  float m_calculated_mother_decaytime_err = -1;
  float m_calculated_mother_decaylength = -1;
  float m_calculated_mother_decaylength_err = -1;
  float m_calculated_mother_dira = -1;
  float m_calculated_mother_fdchi2 = -1;
  float m_calculated_mother_ip = -1;
  float m_calculated_mother_ip_xy = -1;
  float m_calculated_mother_ipchi2 = -1;
  float m_calculated_mother_ip_err = -1;
  float m_calculated_mother_x = -1;
  float m_calculated_mother_y = -1;
  float m_calculated_mother_z = -1;
  float m_calculated_mother_px = -1;
  float m_calculated_mother_py = -1;
  float m_calculated_mother_pz = -1;
  float m_calculated_mother_pe = -1;
  float m_calculated_mother_p = -1;
  float m_calculated_mother_p_err = -1;
  float m_calculated_mother_pt = -1;
  float m_calculated_mother_pt_err = -1;
  int m_calculated_mother_q = -1;
  float m_calculated_mother_eta = -1;
  float m_calculated_mother_rapidity = -1;
  float m_calculated_mother_theta = -1;
  float m_calculated_mother_phi = -1;
  float m_calculated_mother_v = -1;
  float m_calculated_mother_chi2 = -1;
  int m_calculated_mother_ndof = -1;
  int m_calculated_mother_pdgID = -1;
  //float *m_calculated_mother_cov;
  float m_calculated_mother_cov[21] = {0};

  static const int max_intermediates = 8;
  float m_calculated_intermediate_mass[max_intermediates] = {0};
  float m_calculated_intermediate_mass_err[max_intermediates] = {0};
  float m_calculated_intermediate_decaytime[max_intermediates] = {0};
  float m_calculated_intermediate_decaytime_err[max_intermediates] = {0};
  float m_calculated_intermediate_decaylength[max_intermediates] = {0};
  float m_calculated_intermediate_decaylength_err[max_intermediates] = {0};
  float m_calculated_intermediate_dira[max_intermediates] = {0};
  float m_calculated_intermediate_fdchi2[max_intermediates] = {0};
  float m_calculated_intermediate_ip[max_intermediates] = {0};
  float m_calculated_intermediate_ip_xy[max_intermediates] = {0};
  float m_calculated_intermediate_ipchi2[max_intermediates] = {0};
  float m_calculated_intermediate_ip_err[max_intermediates] = {0};
  float m_calculated_intermediate_x[max_intermediates] = {0};
  float m_calculated_intermediate_y[max_intermediates] = {0};
  float m_calculated_intermediate_z[max_intermediates] = {0};
  float m_calculated_intermediate_px[max_intermediates] = {0};
  float m_calculated_intermediate_py[max_intermediates] = {0};
  float m_calculated_intermediate_pz[max_intermediates] = {0};
  float m_calculated_intermediate_pe[max_intermediates] = {0};
  float m_calculated_intermediate_p[max_intermediates] = {0};
  float m_calculated_intermediate_p_err[max_intermediates] = {0};
  float m_calculated_intermediate_pt[max_intermediates] = {0};
  float m_calculated_intermediate_pt_err[max_intermediates] = {0};
  float m_calculated_intermediate_q[max_intermediates] = {0};
  float m_calculated_intermediate_eta[max_intermediates] = {0};
  float m_calculated_intermediate_rapidity[max_intermediates] = {0};
  float m_calculated_intermediate_theta[max_intermediates] = {0};
  float m_calculated_intermediate_phi[max_intermediates] = {0};
  float m_calculated_intermediate_v[max_intermediates] = {0};
  float m_calculated_intermediate_chi2[max_intermediates] = {0};
  float m_calculated_intermediate_ndof[max_intermediates] = {0};
  int m_calculated_intermediate_pdgID[max_intermediates] = {0};
  //float *m_calculated_intermediate_cov[max_intermediates];
  float m_calculated_intermediate_cov[max_intermediates][21] = {{0}, {0}};

  //static const int max_tracks = 20;
  float m_calculated_daughter_mass[max_tracks] = {0};
  float m_calculated_daughter_ip[max_tracks] = {0};
  float m_calculated_daughter_ip_xy[max_tracks] = {0};
  float m_calculated_daughter_ipchi2[max_tracks] = {0};
  float m_calculated_daughter_ip_err[max_tracks] = {0};
  float m_calculated_daughter_x[max_tracks] = {0};
  float m_calculated_daughter_y[max_tracks] = {0};
  float m_calculated_daughter_z[max_tracks] = {0};
  float m_calculated_daughter_px[max_tracks] = {0};
  float m_calculated_daughter_py[max_tracks] = {0};
  float m_calculated_daughter_pz[max_tracks] = {0};
  float m_calculated_daughter_pe[max_tracks] = {0};
  float m_calculated_daughter_p[max_tracks] = {0};
  float m_calculated_daughter_p_err[max_tracks] = {0};
  float m_calculated_daughter_pt[max_tracks] = {0};
  float m_calculated_daughter_pt_err[max_tracks] = {0};
  float m_calculated_daughter_jt[max_tracks] = {0};
  int m_calculated_daughter_q[max_tracks] = {0};
  float m_calculated_daughter_eta[max_tracks] = {0};
  float m_calculated_daughter_rapidity[max_tracks] = {0};
  float m_calculated_daughter_theta[max_tracks] = {0};
  float m_calculated_daughter_phi[max_tracks] = {0};
  float m_calculated_daughter_chi2[max_tracks] = {0};
  int m_calculated_daughter_ndof[max_tracks] = {0};
  int m_calculated_daughter_trid[max_tracks] = {0};
  int m_calculated_daughter_pdgID[max_tracks] = {0};
  //float *m_calculated_daughter_cov[max_tracks];
  float m_calculated_daughter_cov[max_tracks][21] = {{0}, {0}};

  float m_daughter_dca[99] = {0};

  float m_calculated_vertex_x = -1;
  float m_calculated_vertex_y = -1;
  float m_calculated_vertex_z = -1;
  float m_calculated_vertex_v = -1;
  int m_calculated_vertex_nTracks = -1;
  float m_calculated_vertex_chi2 = -1;
  float m_calculated_vertex_ndof = -1;
  //float *m_calculated_vertex_cov;
  float m_calculated_vertex_cov[6] = {0};

  float m_sv_mass = -1;

  int m_nPVs = -1;
  int m_multiplicity = -1;

  int m_runNumber = -1;
  int m_evtNumber = -1;
};

#endif
