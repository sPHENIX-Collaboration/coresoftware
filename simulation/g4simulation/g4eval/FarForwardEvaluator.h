#ifndef G4EVAL_FARFORWARDEVALUATOR_H
#define G4EVAL_FARFORWARDEVALUATOR_H

//===============================================
/// \file FarForwardEvaluator.h
/// \brief Compares reconstructed tracks to truth particles
/// \author Wenliang (Bill) Li
//===============================================

#include <fun4all/SubsysReco.h>

#include <fun4all/Fun4AllHistoManager.h>

#include <set>
#include <string>
#include "TH1.h"
#include "TH2.h"

class CaloEvalStack;
class PHCompositeNode;
class TFile;
class TNtuple;
class TTree;  //Added by Barak

/// \class FarForwardEvaluator
///
/// \brief Compares reconstructed showers to truth particles
///
/// Plan: This module will trace the reconstructed clusters back to
/// the greatest contributor Monte Carlo particle and then
/// test one against the other.
///
class FarForwardEvaluator : public SubsysReco
{
 public:
  FarForwardEvaluator(const std::string &name = "FARFORWARDEVALUATOR",
                const std::string &ffrname = "FFR",
                const std::string &filename = "g4eval_ffr.root",
		const std::string &ip_str = "IP6"
);
  ~FarForwardEvaluator() override{};

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

//
//  void set_strict(bool b) { _strict = b; }
//  // funtions to limit the tracing to only part of the event ---------
//  // and speed up the evaluation
//
//  // when tracing truth showers limit the trace to showers
//  // that result from truth particles with a particular embed flag set
//  // useful if you only want to know about that electron you
//  // embedded into a central hijing event
//  // (evaluation for truth objects with matching embed flag set unaffected)
//  void add_truth_tracing_embed_flag(int flag)
//  {
//    _truth_trace_embed_flags.insert(flag);
//  }
//
//  // limit the tracing of truth particles to those above some
//  // theshold energy. useful for tracing only high energy particles
//  // and ignoring low energy truth particles from a hijing event
//  // (evaluation for objects above threshold unaffected)
//  void set_truth_tracing_energy_threshold(float thresh)
//  {
//    _truth_e_threshold = thresh;
//  }
//
//  // limit the tracing of towers and clusters back to the truth particles
//  // to only those reconstructed objects above a particular energy
//  // threshold (evaluation for objects above threshold unaffected)
//  void set_reco_tracing_energy_threshold(float thresh)
//  {
//    _reco_e_threshold = thresh;
//  }
//
//  // functions to limit the output size ------------------
//  // will no evaluate or write out these particular ntuples
//  // mostly intended for size savings, but some time savings will result
//  void set_do_gpoint_eval(bool b) { _do_gpoint_eval = b; }
//  void set_do_gshower_eval(bool b) { _do_gshower_eval = b; }
//  void set_do_tower_eval(bool b) { _do_tower_eval = b; }
//  void set_do_cluster_eval(bool b) { _do_cluster_eval = b; }
//

private:

  std::string _ffrname;
  std::string _ip_str;

//  TFile *outfile;
//  std::string outfilename;

  TNtuple *g4hitntuple;
  TNtuple *clusterntuple;

  TH2F* h2_ZDC_XY; 
  TH2F* h2_ZDC_XY_double; 

  TH2F* h2_B0_XY; 

  TH2F* h2_RP_XY; 

  Fun4AllHistoManager *hm;

  unsigned int _ievent;

  //Added by Barak
  unsigned int _towerID_debug;
  int _ieta_debug;
  int _iphi_debug;
  float _eta_debug;
  float _phi_debug;
  float _e_debug;
  float _x_debug;
  float _y_debug;
  float _z_debug;

  std::set<int> _truth_trace_embed_flags;
  float _truth_e_threshold;
  float _reco_e_threshold;

//  CaloEvalStack *_caloevalstack;

  //----------------------------------
  // evaluator output ntuples

  bool _strict;

  bool _do_gpoint_eval;
  bool _do_gshower_eval;
  bool _do_tower_eval;
  bool _do_cluster_eval;

  TNtuple *_ntp_gpoint;
  TNtuple *_ntp_gshower;
  TNtuple *_ntp_tower;
  TTree *_tower_debug;  //Added by Barak
  TNtuple *_ntp_cluster;

  int ZDC_hit;
  int event_itt;

  // evaluator output file
  std::string _filename;
  TFile *_tfile;

  // subroutines
  int process_g4hits_ZDC(PHCompositeNode *);
  int process_g4hits_RomanPots(PHCompositeNode *);
  int process_g4hits_B0(PHCompositeNode *);
  
  TH1F* h1_E_dep_smeared;
  TH1F* h1_E_dep;




//
//  // subroutines
//  void printInputInfo(PHCompositeNode *topNode);     ///< print out the input object information (debugging upstream components)
//  void fillOutputNtuples(PHCompositeNode *topNode);  ///< dump the evaluator information into ntuple for external analysis
//  void printOutputInfo(PHCompositeNode *topNode);    ///< print out the ancestry information for detailed diagnosis
//
};

#endif  // G4EVAL_FARFORWARDEVALUATOR_H
