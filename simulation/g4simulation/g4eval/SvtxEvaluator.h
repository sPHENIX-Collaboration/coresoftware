#ifndef G4EVAL_SVTXEVALUATOR_H
#define G4EVAL_SVTXEVALUATOR_H

//===============================================
/// \file SvtxEvaluator.h
/// \brief Compares reconstructed tracks to truth particles
/// \author Michael P. McCumber (revised SVTX version)
//===============================================

#include <fun4all/SubsysReco.h>

#include <string>
#include <set>
#include <vector>

class PHCompositeNode;
class PHTimer;
class SvtxEvalStack;
class TFile;
class TNtuple;
class PHG4Hit;

/// \class SvtxEvaluator
///
/// \brief Compares reconstructed tracks to truth particles
///
/// Plan: This module will trace the reconstructed tracks back to
/// the greatest contributor Monte Carlo particle and then
/// test one against the other.
///
class SvtxEvaluator : public SubsysReco
{
 public:
  SvtxEvaluator(const std::string &name = "SVTXEVALUATOR",
                const std::string &filename = "g4eval.root",
                const std::string &trackmapname = "SvtxTrackMap",
                unsigned int nlayers_maps = 3,
                unsigned int nlayers_intt = 8,
                unsigned int nlayers_tpc = 60);
  virtual ~SvtxEvaluator();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void set_strict(bool b) { _strict = b; }
  void set_use_initial_vertex(bool use_init_vtx) {_use_initial_vertex = use_init_vtx;}
  void do_info_eval(bool b) { _do_info_eval = b; }
  void do_vertex_eval(bool b) { _do_vertex_eval = b; }
  void do_gpoint_eval(bool b) { _do_gpoint_eval = b; }
  void do_g4hit_eval(bool b) { _do_g4hit_eval = b; }
  void do_hit_eval(bool b) { _do_hit_eval = b; }
  void do_cluster_eval(bool b) { _do_cluster_eval = b; }
  void do_g4cluster_eval(bool b) { _do_g4cluster_eval = b; }
  void do_gtrack_eval(bool b) { _do_gtrack_eval = b; }
  void do_track_eval(bool b) { _do_track_eval = b; }
  void do_gseed_eval(bool b) { _do_gseed_eval = b; }

  void do_track_match(bool b) { _do_track_match = b; }
  void do_eval_light(bool b) { _do_eval_light = b; }
  void scan_for_embedded(bool b) { _scan_for_embedded = b; }

 private:
  unsigned int _ievent;
  unsigned int _iseed;

  // eval stack
  SvtxEvalStack *_svtxevalstack;

  //----------------------------------
  // evaluator output ntuples

  bool _strict;
  bool _use_initial_vertex;
  unsigned int _errors;

  bool _do_info_eval;
  bool _do_vertex_eval;
  bool _do_gpoint_eval;
  bool _do_g4hit_eval;
  bool _do_hit_eval;
  bool _do_cluster_eval;
  bool _do_g4cluster_eval;
  bool _do_gtrack_eval;
  bool _do_track_eval;
  bool _do_gseed_eval;

  bool _do_track_match;
  bool _do_eval_light;
  bool _scan_for_embedded;

  unsigned int _nlayers_maps = 3;
  unsigned int _nlayers_intt = 8;
  unsigned int _nlayers_tpc = 60;

  TNtuple *_ntp_info;
  TNtuple *_ntp_vertex;
  TNtuple *_ntp_gpoint;
  TNtuple *_ntp_g4hit;
  TNtuple *_ntp_hit;
  TNtuple *_ntp_cluster;
  TNtuple *_ntp_g4cluster;
  TNtuple *_ntp_gtrack;
  TNtuple *_ntp_track;
  TNtuple *_ntp_gseed;

  // evaluator output file
  std::string _filename;
  //Track map name
  std::string _trackmapname;
  TFile *_tfile;

  PHTimer *_timer;

  //  void LayerClusterG4Particle();

  void G4ClusterSize(PHCompositeNode* topNode, unsigned int layer, std::vector<std::vector<double>> contributing_hits_entry, std::vector<std::vector<double>> contributing_hits_exit, float &g4phisize, float &g4zsize);
  void LayerClusterG4Hits(PHCompositeNode* topNode, std::set<PHG4Hit*> truth_hits, std::vector<PHG4Hit*> &contributing_hits, std::vector<double> &contributing_hits_energy, std::vector<std::vector<double>> &contributing_hits_entry, std::vector<std::vector<double>> &contributing_hits_exit, float layer, float &gx, float &gy, float &gz,  float &gt, float &gedep);
  
  float line_circle_intersection(float x[], float y[], float z[], float radius);

  // output subroutines
  void fillOutputNtuples(PHCompositeNode *topNode);  ///< dump the evaluator information into ntuple for external analysis
  void printInputInfo(PHCompositeNode *topNode);     ///< print out the input object information (debugging upstream components)
  void printOutputInfo(PHCompositeNode *topNode);    ///< print out the ancestry information for detailed diagnosis
};

#endif  // G4EVAL_SVTXEVALUATOR_H
