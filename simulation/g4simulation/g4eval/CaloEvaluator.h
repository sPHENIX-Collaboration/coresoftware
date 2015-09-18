#ifndef __CALOEVALUATOR_H__
#define __CALOEVALUATOR_H__

//===============================================
/// \file CaloEvaluator.h
/// \brief Compares reconstructed tracks to truth particles
/// \author Michael P. McCumber (revised sPHENIX version)
//===============================================

#include <fun4all/SubsysReco.h>
#include <phool/PHCompositeNode.h>

#include <TNtuple.h>
#include <TFile.h>

#include <set>
#include <string>

/// \class CaloEvaluator
///
/// \brief Compares reconstructed showers to truth particles
///
/// Plan: This module will trace the reconstructed clusters back to
/// the greatest contributor Monte Carlo particle and then
/// test one against the other.
///
class CaloEvaluator : public SubsysReco {

 public:
 
  CaloEvaluator(const std::string &name = "CALOEVALUATOR",
		const std::string &caloname = "CEMC",
		const std::string &filename = "g4eval_cemc.root");
  virtual ~CaloEvaluator() {};
		
  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  // forward trace only off of primaries with embed flags set
  void add_truth_tracing_embed_flag(int flag) {
    _truth_trace_embed_flags.insert(flag);
  }

  // backward trace only if reco meets an energy threshold requirement
  void set_reco_tracing_energy_threshold(float thresh) {
    _reco_e_threshold = thresh;
  }

  void set_do_gpoint_eval(bool b) {_do_gpoint_eval = b;}
  void set_do_gshower_eval(bool b) {_do_gshower_eval = b;}
  void set_do_tower_eval(bool b) {_do_tower_eval = b;}
  void set_do_cluster_eval(bool b) {_do_cluster_eval = b;}
  
 private:

  std::string _caloname;
  
  unsigned long _ievent;

  std::set<int> _truth_trace_embed_flags;
  float _reco_e_threshold;

  //----------------------------------
  // evaluator output ntuples

  bool _do_gpoint_eval;
  bool _do_gshower_eval;
  bool _do_tower_eval;
  bool _do_cluster_eval;
  
  TNtuple *_ntp_gpoint;
  TNtuple *_ntp_gshower;
  TNtuple *_ntp_tower;
  TNtuple *_ntp_cluster;

  // evaluator output file
  std::string _filename;
  TFile *_tfile;

  // subroutines
  void printInputInfo(PHCompositeNode *topNode);    ///< print out the input object information (debugging upstream components)
  void fillOutputNtuples(PHCompositeNode *topNode); ///< dump the evaluator information into ntuple for external analysis
  void printOutputInfo(PHCompositeNode *topNode);   ///< print out the ancestry information for detailed diagnosis
};

#endif // __CALOEVALUATOR_H__
