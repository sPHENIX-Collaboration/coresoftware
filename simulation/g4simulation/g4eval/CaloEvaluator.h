#ifndef G4EVAL_CALOEVALUATOR_H
#define G4EVAL_CALOEVALUATOR_H

//===============================================
/// \file CaloEvaluator.h
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
class TTree;  // Added by Barak

/// \class CaloEvaluator
///
/// \brief Compares reconstructed showers to truth particles
///
/// Plan: This module will trace the reconstructed clusters back to
/// the greatest contributor Monte Carlo particle and then
/// test one against the other.
///
class CaloEvaluator : public SubsysReco
{
 public:
  CaloEvaluator(const std::string &name = "CALOEVALUATOR",
                const std::string &caloname = "CEMC",
                const std::string &filename = "g4eval_cemc.root");
  ~CaloEvaluator() override{};

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  // allow user to set the value of the event number
  // useful for condor simulation submissions
  // must be called after Init()
  void set_event(int ievent)
  {
    _ievent = ievent;
  }

  void set_strict(bool b) { _strict = b; }
  // funtions to limit the tracing to only part of the event ---------
  // and speed up the evaluation

  // when tracing truth showers limit the trace to showers
  // that result from truth particles with a particular embed flag set
  // useful if you only want to know about that electron you
  // embedded into a central hijing event
  // (evaluation for truth objects with matching embed flag set unaffected)
  void add_truth_tracing_embed_flag(int flag)
  {
    _truth_trace_embed_flags.insert(flag);
  }

  // limit the tracing of truth particles to those above some
  // theshold energy. useful for tracing only high energy particles
  // and ignoring low energy truth particles from a hijing event
  // (evaluation for objects above threshold unaffected)
  void set_truth_tracing_energy_threshold(float thresh)
  {
    _truth_e_threshold = thresh;
  }

  // limit the tracing of towers and clusters back to the truth particles
  // to only those reconstructed objects above a particular energy
  // threshold (evaluation for objects above threshold unaffected)
  void set_reco_tracing_energy_threshold(float thresh)
  {
    _reco_e_threshold = thresh;
  }

  // functions to limit the output size ------------------
  // will no evaluate or write out these particular ntuples
  // mostly intended for size savings, but some time savings will result
  void set_do_gpoint_eval(bool b) { _do_gpoint_eval = b; }
  void set_do_gshower_eval(bool b) { _do_gshower_eval = b; }
  void set_do_tower_eval(bool b) { _do_tower_eval = b; }
  void set_do_cluster_eval(bool b) { _do_cluster_eval = b; }
  void set_use_towerinfo(bool b) { _use_towerinfo = b; }

 private:
  // subroutines
  void printInputInfo(PHCompositeNode *topNode);     ///< print out the input object information (debugging upstream components)
  void fillOutputNtuples(PHCompositeNode *topNode);  ///< dump the evaluator information into ntuple for external analysis
  void printOutputInfo(PHCompositeNode *topNode);    ///< print out the ancestry information for detailed diagnosis

  CaloEvalStack *_caloevalstack = nullptr;
  TFile *_tfile = nullptr;
  TNtuple *_ntp_cluster = nullptr;
  TNtuple *_ntp_gpoint = nullptr;
  TNtuple *_ntp_gshower = nullptr;
  TNtuple *_ntp_tower = nullptr;
  TTree *_tower_debug = nullptr;  // Added by Barak

  unsigned int _ievent = 0;
  unsigned int _towerID_debug = 0;
  unsigned int m_EvtCounter = 0;

  int _ieta_debug = 0;
  int _iphi_debug = 0;

  float _eta_debug = 0.;
  float _phi_debug = 0.;
  float _e_debug = 0.;
  float _x_debug = 0.;
  float _y_debug = 0.;
  float _z_debug = 0.;
  float _truth_e_threshold = 0.;
  float _reco_e_threshold = 0.;

  bool _do_cluster_eval = true;
  bool _do_gpoint_eval = true;
  bool _do_gshower_eval = true;
  bool _do_tower_eval = true;
  bool _use_towerinfo = false;
  bool _strict = false;

  std::string _caloname;
  std::string _filename;
  std::set<int> _truth_trace_embed_flags;
};

#endif  // G4EVAL_CALOEVALUATOR_H
