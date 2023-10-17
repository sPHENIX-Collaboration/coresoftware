#ifndef G4EVAL_EVENTEVALUATOR_H
#define G4EVAL_EVENTEVALUATOR_H

//===============================================
/// \file EventEvaluator.h
/// \brief Compares reconstructed tracks to truth particles
/// \author Michael P. McCumber (revised sPHENIX version)
//===============================================

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class CaloEvalStack;
class PHCompositeNode;
class PHHepMCGenEventMap;
class PHHepMCGenEvent;
class TFile;
class TNtuple;
class TTree;  // Added by Barak

/// \class EventEvaluator
///
/// \brief Compares reconstructed showers to truth particles
///
/// Plan: This module will trace the reconstructed clusters back to
/// the greatest contributor Monte Carlo particle and then
/// test one against the other.
///
class EventEvaluator : public SubsysReco
{
 public:
  enum class TrackSource_t : unsigned short
  {
    all = 0,
    inner = 1
  };

  EventEvaluator(const std::string& name = "EventEvaluator",
                 const std::string& filename = "g4eval_cemc.root");
  ~EventEvaluator() override{};

  int Init(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  void set_strict(bool b) { _strict = b; }

  void set_do_store_event_level_info(bool b) { _do_store_event_info = b; }
  void set_do_HCALIN(bool b) { _do_HCALIN = b; }
  void set_do_HCALOUT(bool b) { _do_HCALOUT = b; }
  void set_do_CEMC(bool b) { _do_CEMC = b; }
  void set_do_HITS(bool b) { _do_HITS = b; }
  void set_do_TRACKS(bool b) { _do_TRACKS = b; }
  void set_do_CLUSTERS(bool b) { _do_CLUSTERS = b; }
  void set_do_VERTEX(bool b) { _do_VERTEX = b; }
  void set_do_PROJECTIONS(bool b) { _do_PROJECTIONS = b; }
  void set_do_MCPARTICLES(bool b) { _do_MCPARTICLES = b; }
  void set_do_HEPMC(bool b) { _do_HEPMC = b; }
  void set_do_GEOMETRY(bool b) { _do_GEOMETRY = b; }

  // limit the tracing of towers and clusters back to the truth particles
  // to only those reconstructed objects above a particular energy
  // threshold (evaluation for objects above threshold unaffected)
  void set_reco_tracing_energy_threshold(float thresh)
  {
    _reco_e_threshold = thresh;
  }
  void set_reco_tracing_energy_threshold_BECAL(float thresh)
  {
    _reco_e_threshold_BECAL = thresh;
  }

  //! max depth/generation of the MC_particle/PHG4Particle that would be saved.
  void set_depth_MCstack(int d)
  {
    _depth_MCstack = d;
  }

 private:
  bool _do_store_event_info = false;
  bool _do_HCALIN = false;
  bool _do_HCALOUT = false;
  bool _do_CEMC = false;
  bool _do_HITS = false;
  bool _do_TRACKS = false;
  bool _do_CLUSTERS = false;
  bool _do_VERTEX = false;
  bool _do_PROJECTIONS = false;
  bool _do_MCPARTICLES = false;
  bool _do_HEPMC = false;
  bool _do_GEOMETRY = false;
  unsigned int _ievent = 0;

  // Event level info
  float _cross_section = 0.;
  float _event_weight = 0.;
  int _n_generator_accepted = 0;

  // track hits
  int _nHitsLayers = 0;
  int* _hits_layerID = nullptr;
  int* _hits_trueID = nullptr;
  float* _hits_x = nullptr;
  float* _hits_y = nullptr;
  float* _hits_z = nullptr;
  float* _hits_t = nullptr;

  // towers
  int _nTowers_HCALIN = 0;
  float* _tower_HCALIN_E = nullptr;
  int* _tower_HCALIN_iEta = nullptr;
  int* _tower_HCALIN_iPhi = nullptr;
  int* _tower_HCALIN_trueID = nullptr;

  int _nTowers_HCALOUT = 0;
  float* _tower_HCALOUT_E = nullptr;
  int* _tower_HCALOUT_iEta = nullptr;
  int* _tower_HCALOUT_iPhi = nullptr;
  int* _tower_HCALOUT_trueID = nullptr;

  int _nTowers_CEMC = 0;
  float* _tower_CEMC_E = nullptr;
  int* _tower_CEMC_iEta = nullptr;
  int* _tower_CEMC_iPhi = nullptr;
  int* _tower_CEMC_trueID = nullptr;

  // clusters
  int _nclusters_HCALIN = 0;
  float* _cluster_HCALIN_E = nullptr;
  float* _cluster_HCALIN_Eta = nullptr;
  float* _cluster_HCALIN_Phi = nullptr;
  int* _cluster_HCALIN_NTower = nullptr;
  int* _cluster_HCALIN_trueID = nullptr;

  int _nclusters_HCALOUT = 0;
  float* _cluster_HCALOUT_E = nullptr;
  float* _cluster_HCALOUT_Eta = nullptr;
  float* _cluster_HCALOUT_Phi = nullptr;
  int* _cluster_HCALOUT_NTower = nullptr;
  int* _cluster_HCALOUT_trueID = nullptr;

  int _nclusters_CEMC = 0;
  float* _cluster_CEMC_E = nullptr;
  float* _cluster_CEMC_Eta = nullptr;
  float* _cluster_CEMC_Phi = nullptr;
  int* _cluster_CEMC_NTower = nullptr;
  int* _cluster_CEMC_trueID = nullptr;

  // vertex
  float _vertex_x = 0.;
  float _vertex_y = 0.;
  float _vertex_z = 0.;
  int _vertex_NCont = 0;
  float _vertex_true_x = 0.;
  float _vertex_true_y = 0.;
  float _vertex_true_z = 0.;

  // tracks
  int _nTracks = 0;
  float* _track_ID = nullptr;
  float* _track_px = nullptr;
  float* _track_py = nullptr;
  float* _track_pz = nullptr;
  float* _track_dca = nullptr;
  float* _track_dca_2d = nullptr;
  float* _track_trueID = nullptr;
  unsigned short* _track_source = nullptr;

  int _nProjections = 0;
  float* _track_ProjTrackID = nullptr;
  int* _track_ProjLayer = nullptr;
  float* _track_TLP_x = nullptr;
  float* _track_TLP_y = nullptr;
  float* _track_TLP_z = nullptr;
  float* _track_TLP_t = nullptr;
  float* _track_TLP_true_x = nullptr;
  float* _track_TLP_true_y = nullptr;
  float* _track_TLP_true_z = nullptr;
  float* _track_TLP_true_t = nullptr;

  // MC particles
  int _nMCPart = 0;
  int* _mcpart_ID = nullptr;
  int* _mcpart_ID_parent = nullptr;
  int* _mcpart_PDG = nullptr;
  float* _mcpart_E = nullptr;
  float* _mcpart_px = nullptr;
  float* _mcpart_py = nullptr;
  float* _mcpart_pz = nullptr;
  int* _mcpart_BCID = nullptr;

  // MC particles
  int _nHepmcp = 0;
  int _hepmcp_procid = 0;
  float _hepmcp_x1 = 0.;
  float _hepmcp_x2 = 0.;
  //  float* _hepmcp_ID_parent;
  int* _hepmcp_status = nullptr;
  int* _hepmcp_PDG = nullptr;
  float* _hepmcp_E = nullptr;
  float* _hepmcp_px = nullptr;
  float* _hepmcp_py = nullptr;
  float* _hepmcp_pz = nullptr;
  int* _hepmcp_m1 = nullptr;
  int* _hepmcp_m2 = nullptr;
  int* _hepmcp_BCID = nullptr;

  int _calo_ID = 0;
  int _calo_towers_N = 0;
  int* _calo_towers_iEta = nullptr;
  int* _calo_towers_iPhi = nullptr;
  float* _calo_towers_Eta = nullptr;
  float* _calo_towers_Phi = nullptr;
  float* _calo_towers_x = nullptr;
  float* _calo_towers_y = nullptr;
  float* _calo_towers_z = nullptr;
  int* _geometry_done = nullptr;

  float _reco_e_threshold = 0.;
  float _reco_e_threshold_BECAL = 0.;
  int _depth_MCstack = 0;

  CaloEvalStack* _caloevalstackHCALIN = nullptr;
  CaloEvalStack* _caloevalstackHCALOUT = nullptr;
  CaloEvalStack* _caloevalstackCEMC = nullptr;

  //----------------------------------
  // evaluator output ntuples

  bool _strict = false;

  TTree* _event_tree = nullptr;     // Added by Barak
  TTree* _geometry_tree = nullptr;  // Added by Barak

  // evaluator output file
  std::string _filename;
  TFile* _tfile = nullptr;
  TFile* _tfile_geometry = nullptr;

  // subroutines
  int GetProjectionIndex(const std::string& projname);    ///< return track projection index for given track projection layer
  std::string GetProjectionNameFromIndex(int projindex);  ///< return track projection layer name from projection index (see GetProjectionIndex)
  void fillOutputNtuples(PHCompositeNode* topNode);       ///< dump the evaluator information into ntuple for external analysis
  void resetGeometryArrays();                             ///< reset the tree variables before filling for a new event
  void resetBuffer();                                     ///< reset the tree variables before filling for a new event

  const int _maxNHits = 10000;
  const int _maxNTowersCentral = 2000;
  const int _maxNTowersCalo = 5000000;
  const int _maxNclustersCentral = 2000;
  const int _maxNTracks = 200;
  const int _maxNProjections = 2000;
  const int _maxNMCPart = 100000;
  const int _maxNHepmcp = 1000;

  enum calotype
  {
    kCEMC = 0,
    kHCALIN = 1,
    kHCALOUT = 2
  };
};

#endif  // G4EVAL_EVENTEVALUATOR_H
