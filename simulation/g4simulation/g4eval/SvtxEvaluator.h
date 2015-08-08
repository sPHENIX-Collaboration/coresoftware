#ifndef __SVTXEVALUATOR_H__
#define __SVTXEVALUATOR_H__

//===============================================
/// \file SvtxEvaluator.h
/// \brief Compares reconstructed tracks to truth particles
/// \author Michael P. McCumber (revised SVTX version)
//===============================================

// PHENIX includes
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>

#include <TNtuple.h>
#include <TFile.h>

#include <string>

/// \class SvtxEvaluator
///
/// \brief Compares reconstructed tracks to truth particles
///
/// Plan: This module will trace the reconstructed tracks back to
/// the greatest contributor Monte Carlo particle and then
/// test one against the other.
///
class SvtxEvaluator : public SubsysReco {
  
public:
 
  SvtxEvaluator(const std::string &name = "SVTXEVALUATOR",
                const std::string &filename = "g4eval.root");
  virtual ~SvtxEvaluator() {};
		
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void do_vertex_eval(bool b) {_do_vertex_eval = b;}
  void do_gpoint_eval(bool b) {_do_gpoint_eval = b;}
  void do_g4hit_eval(bool b) {_do_g4hit_eval = b;}
  void do_hit_eval(bool b) {_do_hit_eval = b;}
  void do_cluster_eval(bool b) {_do_cluster_eval = b;}
  void do_gtrack_eval(bool b) {_do_gtrack_eval = b;}
  void do_track_eval(bool b) {_do_track_eval = b;}
  
 private:

  unsigned int _ievent;

  //----------------------------------
  // evaluator output ntuples

  bool _do_vertex_eval;
  bool _do_gpoint_eval;
  bool _do_g4hit_eval;
  bool _do_hit_eval;
  bool _do_cluster_eval;
  bool _do_gtrack_eval;
  bool _do_track_eval;
  
  TNtuple *_ntp_vertex;
  TNtuple *_ntp_gpoint;
  TNtuple *_ntp_g4hit;
  TNtuple *_ntp_hit;
  TNtuple *_ntp_cluster;
  TNtuple *_ntp_gtrack;
  TNtuple *_ntp_track;

  // evaluator output file
  std::string _filename;
  TFile *_tfile;

  // output subroutines
  void fillOutputNtuples(PHCompositeNode* topNode); ///< dump the evaluator information into ntuple for external analysis
  void printInputInfo(PHCompositeNode* topNode);    ///< print out the input object information (debugging upstream components)
  void printOutputInfo(PHCompositeNode* topNode);   ///< print out the ancestry information for detailed diagnosis
};

#endif // __SVTXEVALUATOR_H__
