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
class TFile;
class TNtuple;
class TTree;  //Added by Barak

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
  EventEvaluator(const std::string& name = "EventEvaluator",
                 const std::string& filename = "g4eval_cemc.root");
  virtual ~EventEvaluator(){};

  int Init(PHCompositeNode* topNode);
  int process_event(PHCompositeNode* topNode);
  int End(PHCompositeNode* topNode);

  void set_strict(bool b) { _strict = b; }

  void set_do_FHCAL(bool b) { _do_FHCAL = b; }
  void set_do_HCALIN(bool b) { _do_HCALIN = b; }
  void set_do_HCALOUT(bool b) { _do_HCALOUT = b; }
  void set_do_EHCAL(bool b) { _do_EHCAL = b; }
  void set_do_FEMC(bool b) { _do_FEMC = b; }
  void set_do_CEMC(bool b) { _do_CEMC = b; }
  void set_do_EEMC(bool b) { _do_EEMC = b; }
  void set_do_DRCALO(bool b) { _do_DRCALO = b; }
  void set_do_HITS(bool b) { _do_HITS = b; }
  void set_do_TRACKS(bool b) { _do_TRACKS = b; }
  void set_do_CLUSTERS(bool b) { _do_CLUSTERS = b; }
  void set_do_VERTEX(bool b) { _do_VERTEX = b; }
  void set_do_PROJECTIONS(bool b) { _do_PROJECTIONS = b; }
  void set_do_MCPARTICLES(bool b) { _do_MCPARTICLES = b; }
  // funtions to limit the tracing to only part of the event ---------
  // and speed up the evaluation

  // limit the tracing of towers and clusters back to the truth particles
  // to only those reconstructed objects above a particular energy
  // threshold (evaluation for objects above threshold unaffected)
  void set_reco_tracing_energy_threshold(float thresh)
  {
    _reco_e_threshold = thresh;
  }

private:
  bool _do_FHCAL;
  bool _do_HCALIN;
  bool _do_HCALOUT;
  bool _do_EHCAL;
  bool _do_FEMC;
  bool _do_CEMC;
  bool _do_EEMC;
  bool _do_DRCALO;
  bool _do_HITS;
  bool _do_TRACKS;
  bool _do_CLUSTERS;
  bool _do_VERTEX;
  bool _do_PROJECTIONS;
  bool _do_MCPARTICLES;
  unsigned int _ievent;

  // track hits
  int _nHitsLayers;
  int* _hits_layerID;
  int* _hits_trueID;
  float* _hits_x;
  float* _hits_y;
  float* _hits_z;
  float* _hits_t;

  // towers
  int _nTowers_FHCAL;
  float* _tower_FHCAL_E;
  int* _tower_FHCAL_iEta;
  int* _tower_FHCAL_iPhi;
  int* _tower_FHCAL_trueID;

  // towers
  int _nTowers_HCALIN;
  float* _tower_HCALIN_E;
  int* _tower_HCALIN_iEta;
  int* _tower_HCALIN_iPhi;
  int* _tower_HCALIN_trueID;

  // towers
  int _nTowers_HCALOUT;
  float* _tower_HCALOUT_E;
  int* _tower_HCALOUT_iEta;
  int* _tower_HCALOUT_iPhi;
  int* _tower_HCALOUT_trueID;

  int _nTowers_EHCAL;
  float* _tower_EHCAL_E;
  int* _tower_EHCAL_iEta;
  int* _tower_EHCAL_iPhi;
  int* _tower_EHCAL_trueID;
  
  int _nTowers_DRCALO;
  float* _tower_DRCALO_E;
  int* _tower_DRCALO_NScint;
  int* _tower_DRCALO_NCerenkov;
  int* _tower_DRCALO_iEta;
  int* _tower_DRCALO_iPhi;
  int* _tower_DRCALO_trueID;

  int _nTowers_FEMC;
  float* _tower_FEMC_E;
  int* _tower_FEMC_iEta;
  int* _tower_FEMC_iPhi;
  int* _tower_FEMC_trueID;

  int _nTowers_EEMC;
  float* _tower_EEMC_E;
  int* _tower_EEMC_iEta;
  int* _tower_EEMC_iPhi;
  int* _tower_EEMC_trueID;

  int _nTowers_CEMC;
  float* _tower_CEMC_E;
  int* _tower_CEMC_iEta;
  int* _tower_CEMC_iPhi;
  int* _tower_CEMC_trueID;
  
  // clusters
  int _nclusters_FHCAL;
  float* _cluster_FHCAL_E;
  float* _cluster_FHCAL_Eta;
  float* _cluster_FHCAL_Phi;
  int* _cluster_FHCAL_NTower;
  int* _cluster_FHCAL_trueID;

  int _nclusters_HCALIN;
  float* _cluster_HCALIN_E;
  float* _cluster_HCALIN_Eta;
  float* _cluster_HCALIN_Phi;
  int* _cluster_HCALIN_NTower;
  int* _cluster_HCALIN_trueID;

  int _nclusters_HCALOUT;
  float* _cluster_HCALOUT_E;
  float* _cluster_HCALOUT_Eta;
  float* _cluster_HCALOUT_Phi;
  int* _cluster_HCALOUT_NTower;
  int* _cluster_HCALOUT_trueID;

  int _nclusters_EHCAL;
  float* _cluster_EHCAL_E;
  float* _cluster_EHCAL_Eta;
  float* _cluster_EHCAL_Phi;
  int* _cluster_EHCAL_NTower;
  int* _cluster_EHCAL_trueID;
  
  int _nclusters_FEMC;
  float* _cluster_FEMC_E;
  float* _cluster_FEMC_Eta;
  float* _cluster_FEMC_Phi;
  int* _cluster_FEMC_NTower;
  int* _cluster_FEMC_trueID;

  int _nclusters_CEMC;
  float* _cluster_CEMC_E;
  float* _cluster_CEMC_Eta;
  float* _cluster_CEMC_Phi;
  int* _cluster_CEMC_NTower;
  int* _cluster_CEMC_trueID;

  int _nclusters_EEMC;
  float* _cluster_EEMC_E;
  float* _cluster_EEMC_Eta;
  float* _cluster_EEMC_Phi;
  int* _cluster_EEMC_NTower;
  int* _cluster_EEMC_trueID;
  
  // vertex
  float _vertex_x;
  float _vertex_y;
  float _vertex_z;
  int _vertex_NCont;
  float _vertex_true_x;
  float _vertex_true_y;
  float _vertex_true_z;

  // tracks
  int _nTracks;
  float* _track_ID;
  float* _track_px;
  float* _track_py;
  float* _track_pz;
  float* _track_trueID;

  int _nProjections;
  float* _track_ProjTrackID;
  int* _track_ProjLayer;
  float* _track_TLP_x;
  float* _track_TLP_y;
  float* _track_TLP_z;
  float* _track_TLP_t;
  float* _track_TLP_true_x;
  float* _track_TLP_true_y;
  float* _track_TLP_true_z;
  float* _track_TLP_true_t;

  // MC particles
  int _nMCPart;
  float* _mcpart_ID;
  float* _mcpart_ID_parent;
  float* _mcpart_PDG;
  float* _mcpart_E;
  float* _mcpart_px;
  float* _mcpart_py;
  float* _mcpart_pz;

  float _reco_e_threshold;
  int _depth_MCstack;

  CaloEvalStack *_caloevalstackFHCAL;
  CaloEvalStack *_caloevalstackHCALIN;
  CaloEvalStack *_caloevalstackHCALOUT;
  CaloEvalStack *_caloevalstackEHCAL;
  CaloEvalStack *_caloevalstackDRCALO;
  CaloEvalStack *_caloevalstackFEMC;
  CaloEvalStack *_caloevalstackCEMC;
  CaloEvalStack *_caloevalstackEEMC;

  //----------------------------------
  // evaluator output ntuples

  bool _strict;

  TTree* _event_tree;  //Added by Barak

  // evaluator output file
  std::string _filename;
  TFile* _tfile;

  // subroutines
  int GetProjectionIndex(std::string projname);           ///< return track projection index for given track projection layer
  std::string GetProjectionNameFromIndex(int projindex);  ///< return track projection layer name from projection index (see GetProjectionIndex)
  void fillOutputNtuples(PHCompositeNode* topNode);       ///< dump the evaluator information into ntuple for external analysis
  void resetBuffer();                                     ///< reset the tree variables before filling for a new event

  const int _maxNHits = 5000;
  const int _maxNTowers = 50 * 50;
  const int _maxNTowersDR = 3000 * 3000;
  const int _maxNclusters = 100;
  const int _maxNTracks = 200;
  const int _maxNProjections = 2000;
  const int _maxNMCPart = 10000;
};

#endif  // G4EVAL_EVENTEVALUATOR_H
