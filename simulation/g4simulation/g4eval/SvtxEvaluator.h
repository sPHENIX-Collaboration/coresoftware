#ifndef SVTXEVALUATOR_H__
#define SVTXEVALUATOR_H__

//===============================================
/// \file SvtxEvaluator.h
/// \brief Compares reconstructed tracks to truth particles
/// \author Michael P. McCumber (revised SVTX version)
//===============================================


#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <phool/PHTimer.h>
#include <string>

class PHCompositeNode;

class SvtxEvalStack;
class TFile;
class TNtuple;

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
                const std::string &filename = "g4eval.root",
                const std::string &trackmapname = "SvtxTrackMap",
		unsigned int nlayers_maps = 3,
		unsigned int nlayers_intt = 8,
		unsigned int nlayers_tpc = 60);
  virtual ~SvtxEvaluator() {}
		
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void set_strict(bool b) {_strict = b;}
  
  void do_vertex_eval(bool b) {_do_vertex_eval = b;}
  void do_gpoint_eval(bool b) {_do_gpoint_eval = b;}
  void do_g4hit_eval(bool b) {_do_g4hit_eval = b;}
  void do_hit_eval(bool b) {_do_hit_eval = b;}
  void do_cluster_eval(bool b) {_do_cluster_eval = b;}
  void do_gtrack_eval(bool b) {_do_gtrack_eval = b;}
  void do_track_eval(bool b) {_do_track_eval = b;}
  void do_gseed_eval(bool b) {_do_gseed_eval = b;}

  void do_track_match(bool b) {_do_track_match = b;}
  void do_eval_light(bool b) {_do_eval_light = b;}
  void scan_for_embedded(bool b) {_scan_for_embedded = b;}
  
 private:

  unsigned int _ievent;

  // eval stack
  SvtxEvalStack* _svtxevalstack;
  
  //----------------------------------
  // evaluator output ntuples

  bool _strict;
  unsigned int _errors;
  
  bool _do_vertex_eval;
  bool _do_gpoint_eval;
  bool _do_g4hit_eval;
  bool _do_hit_eval;
  bool _do_cluster_eval;
  bool _do_gtrack_eval;
  bool _do_track_eval;
  bool _do_gseed_eval;

  bool _do_track_match;
  bool _do_eval_light;
  bool _scan_for_embedded;

  unsigned int _nlayers_maps = 3;
  unsigned int _nlayers_intt = 8;
  unsigned int _nlayers_tpc = 60;

  TNtuple *_ntp_vertex;
  TNtuple *_ntp_gpoint;
  TNtuple *_ntp_g4hit;
  TNtuple *_ntp_hit;
  TNtuple *_ntp_cluster;
  TNtuple *_ntp_gtrack;
  TNtuple *_ntp_track;
  TNtuple *_ntp_gseed;

  // evaluator output file
  std::string _filename;
  //Track map name
  std::string _trackmapname;
  TFile *_tfile;

  PHTimer *_timer;

  float line_circle_intersection(float x[], float y[], float z[], float radius);

  // output subroutines
  void fillOutputNtuples(PHCompositeNode* topNode); ///< dump the evaluator information into ntuple for external analysis
  void printInputInfo(PHCompositeNode* topNode);    ///< print out the input object information (debugging upstream components)
  void printOutputInfo(PHCompositeNode* topNode);   ///< print out the ancestry information for detailed diagnosis
};

#endif // __SVTXEVALUATOR_H__
