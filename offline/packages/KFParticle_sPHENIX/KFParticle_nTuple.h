#ifndef KFPARTICLESPHENIX_KFPARTICLENTUPLE_H
#define KFPARTICLESPHENIX_KFPARTICLENTUPLE_H

#include "KFParticle_Config.h"
#include "KFParticle_caloInfo.h"
#include "KFParticle_detectorInfo.h"
#include "KFParticle_triggerInfo.h"
#include "KFParticle_truthInfo.h"
#include "KFParticle_utilities.h"

#include <KFParticle.h>

#include <limits>
#include <string>  // for string
#include <vector>
#include <cstdint>  // fixed width integer types used for BCO counters

class PHCompositeNode;
class TTree;
class CentralityInfo;
class KFParticle_Tools;

class KFParticle_nTuple
{
 public:
  /// Constructor; binds the configuration only (the owner's instance is not constructed yet)
  explicit KFParticle_nTuple(const KFParticle_Config &config)
    : m_config(config)
    , m_detector(config)
    , m_truth(config)
    , m_calo(config)
  {
  }

  /// Destructor
  virtual ~KFParticle_nTuple() = default;

  /// Unused for now, variables are initialised in the header
  void initializeVariables();

  /// Initialises required branches based off the user selection (number of tracks, PV constraints etc ) and sets branch names if specified
  void initializeBranches(PHCompositeNode *topNode);

  /// Fills required information for your selection, also requests truth and detector information if needed
  void fillBranch(PHCompositeNode *topNode,
                  KFParticle motherParticle,
                  const KFParticle &vertex,
                  std::vector<KFParticle> daughters,
                  std::vector<KFParticle> intermediates);

  // pass event-level BCO values from KFParticle_sPHENIX
  void set_event_bcos(const int64_t this_bco, const int64_t last_bco)
  {
    m_event_bco = this_bco;
    m_last_event_bco = last_bco;
  }

  float calc_secondary_vertex_mass_noPID(std::vector<KFParticle> kfp_daughters) const;

  bool fillConditionMet() const;

  /// Set by KFParticle_sPHENIX to itself, for the calculations and node caches that still live in KFParticle_Tools
  void setToolsInstance(KFParticle_Tools *tools) { m_tools = tools; }

 protected:
  bool m_extrapolateTracksToSV_nTuple{true};
  bool isTrackEMCalmatch{true};
  bool m_use_fake_pv_nTuple{false};
  std::string m_qa_tree_path{"DecayTree"};  // name (optionally dir/name) of the tree in the QA histogram manager
  int m_verbosity_nTuple{0};

 private:
  static constexpr int max_tracks{KFParticle_utilities::max_tracks};

  const KFParticle_Config &m_config;

  // The KFParticle_sPHENIX object itself (see setToolsInstance). Non-const because makeAllPrimaryVertices
  // and get_dEdx update its per-event node pointers.
  KFParticle_Tools *m_tools{nullptr};

  KFParticle_detectorInfo m_detector;
  KFParticle_truthInfo m_truth;
  KFParticle_caloInfo m_calo;
  KFParticle_triggerInfo m_trigger;

  // This fills a reduced output meant for the ML momentum calibration
  // It contains mother mass/pT/eta/phi/PID, daughter pT/pTErr/eta/phi/charge/PID/covariance and runNumber only
  void initializeMomentumCalibrationBranches(const std::string &mother_name);
  void fillMomentumCalibrationBranch(PHCompositeNode *topNode,
                                     KFParticle &motherParticle,
                                     KFParticle *daughterArray);

  // Distance of the mother, daughters and intermediates to every primary vertex in the event
  void allPVInfo(PHCompositeNode *topNode,
                 const KFParticle &motherParticle,
                 std::vector<KFParticle> daughters,
                 std::vector<KFParticle> intermediates);
  void clearAllPVInfo();

  // Event-level counts and jT, moved from KFParticle_Tools (only the nTuple uses them)
  void getTracksFromBC(PHCompositeNode *topNode, const int &bunch_crossing, const std::string &vertexMapName, int &nTracks, int &nPVs) const;
  int getTracksFromVertex(PHCompositeNode *topNode, const KFParticle &vertex, const std::string &vertexMapName) const;
  static int getNchargedSiSeedMultiplicity(PHCompositeNode *topNode, const int &bunch_crossing, int verbosity);
  static float calculateJT(const KFParticle &mother, const KFParticle &daughter);


  TTree *m_tree{nullptr};

  float m_calculated_mother_mass{-1};
  float m_calculated_mother_mass_err{-1};
  float m_calculated_mother_decaytime{-1};
  float m_calculated_mother_decaytime_err{-1};
  float m_calculated_mother_decaylength{-1};
  float m_calculated_mother_decaylength_err{-1};
  float m_calculated_mother_decaylength_xy{-1};
  float m_calculated_mother_decaylength_xy_err{-1};
  float m_calculated_mother_dira{-1};
  float m_calculated_mother_dira_xy{-1};
  float m_calculated_mother_fdchi2{-1};
  float m_calculated_mother_PV_dca{-1};
  float m_calculated_mother_PV_dca_xy{-1};
  float m_calculated_mother_PV_dca_sig{-1};
  float m_calculated_mother_PV_dca_err{-1};
  float m_calculated_mother_x{-1};
  float m_calculated_mother_y{-1};
  float m_calculated_mother_z{-1};
  float m_calculated_mother_px{-1};
  float m_calculated_mother_py{-1};
  float m_calculated_mother_pz{-1};
  float m_calculated_mother_pe{-1};
  float m_calculated_mother_p{-1};
  float m_calculated_mother_p_err{-1};
  float m_calculated_mother_pt{-1};
  float m_calculated_mother_pt_err{-1};
  char m_calculated_mother_q{-1};
  float m_calculated_mother_eta{-1};
  float m_calculated_mother_rapidity{-1};
  float m_calculated_mother_theta{-1};
  float m_calculated_mother_phi{-1};
  float m_calculated_mother_v{-1};
  float m_calculated_mother_chi2{-1};
  int m_calculated_mother_ndof{-1};
  int m_calculated_mother_pdgID{-1};
  // float *m_calculated_mother_cov;
  float m_calculated_mother_cov[21]{0};

  static const int max_intermediates{8};
  float m_calculated_intermediate_mass[max_intermediates]{0};
  float m_calculated_intermediate_mass_err[max_intermediates]{0};
  float m_calculated_intermediate_decaytime[max_intermediates]{0};
  float m_calculated_intermediate_decaytime_err[max_intermediates]{0};
  float m_calculated_intermediate_decaylength[max_intermediates]{0};
  float m_calculated_intermediate_decaylength_err[max_intermediates]{0};
  float m_calculated_intermediate_decaylength_xy[max_intermediates]{0};
  float m_calculated_intermediate_decaylength_xy_err[max_intermediates]{0};
  float m_calculated_intermediate_dira[max_intermediates]{0};
  float m_calculated_intermediate_fdchi2[max_intermediates]{0};
  float m_calculated_intermediate_PV_dca[max_intermediates]{0};
  float m_calculated_intermediate_PV_dca_xy[max_intermediates]{0};
  float m_calculated_intermediate_PV_dca_sig[max_intermediates]{0};
  float m_calculated_intermediate_PV_dca_err[max_intermediates]{0};
  float m_calculated_intermediate_x[max_intermediates]{0};
  float m_calculated_intermediate_y[max_intermediates]{0};
  float m_calculated_intermediate_z[max_intermediates]{0};
  float m_calculated_intermediate_px[max_intermediates]{0};
  float m_calculated_intermediate_py[max_intermediates]{0};
  float m_calculated_intermediate_pz[max_intermediates]{0};
  float m_calculated_intermediate_pe[max_intermediates]{0};
  float m_calculated_intermediate_p[max_intermediates]{0};
  float m_calculated_intermediate_p_err[max_intermediates]{0};
  float m_calculated_intermediate_pt[max_intermediates]{0};
  float m_calculated_intermediate_pt_err[max_intermediates]{0};
  char m_calculated_intermediate_q[max_intermediates]{0};
  float m_calculated_intermediate_eta[max_intermediates]{0};
  float m_calculated_intermediate_rapidity[max_intermediates]{0};
  float m_calculated_intermediate_theta[max_intermediates]{0};
  float m_calculated_intermediate_phi[max_intermediates]{0};
  float m_calculated_intermediate_v[max_intermediates]{0};
  float m_calculated_intermediate_chi2[max_intermediates]{0};
  int m_calculated_intermediate_ndof[max_intermediates]{0};
  int m_calculated_intermediate_pdgID[max_intermediates]{0};
  // float *m_calculated_intermediate_cov[max_intermediates];
  float m_calculated_intermediate_cov[max_intermediates][21]{{0}, {0}};

  float m_calculated_daughter_mass[max_tracks]{0};
  float m_calculated_daughter_PV_dca[max_tracks]{0};
  float m_calculated_daughter_PV_dca_err[max_tracks]{0};
  float m_calculated_daughter_PV_dca_sig[max_tracks]{0};
  float m_calculated_daughter_PV_dca_xy[max_tracks]{0};
  float m_calculated_daughter_PV_dca_xy_sig[max_tracks]{0};
  float m_calculated_daughter_x[max_tracks]{0};
  float m_calculated_daughter_y[max_tracks]{0};
  float m_calculated_daughter_z[max_tracks]{0};
  float m_calculated_daughter_px[max_tracks]{0};
  float m_calculated_daughter_py[max_tracks]{0};
  float m_calculated_daughter_pz[max_tracks]{0};
  float m_calculated_daughter_pe[max_tracks]{0};
  float m_calculated_daughter_p[max_tracks]{0};
  float m_calculated_daughter_p_err[max_tracks]{0};
  float m_calculated_daughter_pt[max_tracks]{0};
  float m_calculated_daughter_pt_err[max_tracks]{0};
  float m_calculated_daughter_jt[max_tracks]{0};
  char m_calculated_daughter_q[max_tracks]{0};
  int m_calculated_daughter_bunch_crossing[max_tracks]{0};
  float m_calculated_daughter_eta[max_tracks]{0};
  float m_calculated_daughter_rapidity[max_tracks]{0};
  float m_calculated_daughter_theta[max_tracks]{0};
  float m_calculated_daughter_phi[max_tracks]{0};
  float m_calculated_daughter_chi2[max_tracks]{0};
  int m_calculated_daughter_ndof[max_tracks]{0};
  int m_calculated_daughter_trid[max_tracks]{0};
  int m_calculated_daughter_pdgID[max_tracks]{0};
  // float *m_calculated_daughter_cov[max_tracks];
  float m_calculated_daughter_cov[max_tracks][21]{{0}, {0}};
  float m_calculated_daughter_dedx[max_tracks]{0};
  // float m_calculated_daughter_expected_dedx_pion[max_tracks] {0};
  // float m_calculated_daughter_expected_dedx_kaon[max_tracks] {0};
  // float m_calculated_daughter_expected_dedx_proton[max_tracks] {0}; 
  float m_daughter_dca[99]{0};
  float m_daughter_dca_xy[99]{0};
  float m_daughter_dca_sig[99]{0};
  float m_daughter_dca_sig_xy[99]{0};

  float m_calculated_vertex_x{-1};
  float m_calculated_vertex_y{-1};
  float m_calculated_vertex_z{-1};
  float m_calculated_vertex_v{-1};
  float m_calculated_vertex_chi2{-1};
  unsigned int m_calculated_vertex_ndof{std::numeric_limits<unsigned int>::max()};
  int m_calculated_vertex_ID{-1};
  // float *m_calculated_vertex_cov;
  float m_calculated_vertex_cov[6]{0};
  float m_centrality_mbd{-1};
  CentralityInfo *m_CentInfo{nullptr};

  float m_sv_mass{-1};

  int m_nPVs{-1};
  int m_multiplicity{-1};
  int m_ncharged_siseed_multiplicity{-1};
  int m_nTracksOfVertex{-1};

  int m_runNumber{-1};
  int m_evtNumber{-1};
  int64_t m_bco{-1};
  uint64_t m_event_bco{0};//current event BCO
  uint64_t m_last_event_bco{0}; //only keeping this, BCO for the last event

  std::vector<float> allPV_x;
  std::vector<float> allPV_y;
  std::vector<float> allPV_z;
  std::vector<float> allPV_mother_PV_DCA;
  std::vector<float> allPV_mother_PV_DCA_StdDev;
  std::vector<float> allPV_daughter_PV_DCA[max_tracks];
  std::vector<float> allPV_daughter_PV_DCA_StdDev[max_tracks];
  std::vector<float> allPV_intermediates_PV_DCA[max_tracks];
  std::vector<float> allPV_intermediates_PV_DCA_StdDev[max_tracks];
};

#endif
