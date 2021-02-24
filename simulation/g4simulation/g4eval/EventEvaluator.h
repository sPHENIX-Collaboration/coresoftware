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

#include <TLorentzVector.h>  // for TLorentzVector

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
  EventEvaluator(const std::string &name = "EventEvaluator",
                const std::string &filename = "g4eval_cemc.root");
  virtual ~EventEvaluator(){};

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void set_strict(bool b) { _strict = b; }
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

  unsigned int _ievent;

  // track hits
  int _nHitsLayers;
  TLorentzVector* _hits_TLH_xyzt;
  int* _hits_layerID;

  // towers
  int _nTowers_FHCAL;
  float* _tower_FHCAL_E;
  int* _tower_FHCAL_iEta;
  int* _tower_FHCAL_iPhi;
  float* _tower_FHCAL_trueID;

  // towers
  int _nTowers_FEMC;
  float* _tower_FEMC_E;
  int* _tower_FEMC_iEta;
  int* _tower_FEMC_iPhi;
  float* _tower_FEMC_trueID;

  // vertex
  int _vertex_x;
  int _vertex_y;
  int _vertex_z;

  // tracks
  int _nTracks;
  float* _track_ID;
  float* _track_px;
  float* _track_py;
  float* _track_pz;
  float* _track_trueID;
  TLorentzVector* _track_TLP_FTTL_0;
  TLorentzVector* _track_TLP_FTTL_0_true;
  TLorentzVector* _track_TLP_FTTL_1;
  TLorentzVector* _track_TLP_FTTL_1_true;
  TLorentzVector* _track_TLP_FTTL_2;
  TLorentzVector* _track_TLP_FTTL_2_true;
  TLorentzVector* _track_TLP_ETTL_0;
  TLorentzVector* _track_TLP_ETTL_0_true;
  TLorentzVector* _track_TLP_ETTL_1;
  TLorentzVector* _track_TLP_ETTL_1_true;
  TLorentzVector* _track_TLP_FHCAL_0;
  TLorentzVector* _track_TLP_FHCAL_0_true;
  TLorentzVector* _track_TLP_FEMC_0;
  TLorentzVector* _track_TLP_FEMC_0_true;
  TLorentzVector* _track_TLP_CTTL_0;
  TLorentzVector* _track_TLP_CTTL_0_true;
  TLorentzVector* _track_TLP_CTTL_1;
  TLorentzVector* _track_TLP_CTTL_1_true;
  TLorentzVector* _track_TLP_CTTL_2;
  TLorentzVector* _track_TLP_CTTL_2_true;

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

  CaloEvalStack *_caloevalstack;

  //----------------------------------
  // evaluator output ntuples

  bool _strict;

  TNtuple *_ntp_gpoint;
  // TNtuple *_ntp_gshower;
  // TNtuple *_ntp_tower;
  TTree *_event_tree;  //Added by Barak
  // TNtuple *_ntp_cluster;

  // evaluator output file
  std::string _filename;
  TFile *_tfile;

  const int _maxNHits = 5000;
  const int _maxNTowers = 50*50;
  const int _maxNProjections = 10;
  const int _maxNTracks = 200;
  const int _maxNMCPart = 1000;

  // subroutines
  int GetProjectionIndex(std::string projname);     ///< print out the input object information (debugging upstream components)
  std::string GetProjectionNameFromIndex(int projindex);
  void printInputInfo(PHCompositeNode *topNode);     ///< print out the input object information (debugging upstream components)
  void fillOutputNtuples(PHCompositeNode *topNode);  ///< dump the evaluator information into ntuple for external analysis
  void printOutputInfo(PHCompositeNode *topNode);    ///< print out the ancestry information for detailed diagnosis
  void resetBuffer();    ///< print out the ancestry information for detailed diagnosis
  void FillTrackProjVector(int trackindex, int trackStateIndex, bool truevalues,float xpostrk,float ypostrk,float zpostrk,float timetrk);

};

#endif  // G4EVAL_EVENTEVALUATOR_H
