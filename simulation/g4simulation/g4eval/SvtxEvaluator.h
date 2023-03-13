#ifndef G4EVAL_SVTXEVALUATOR_H
#define G4EVAL_SVTXEVALUATOR_H

//===============================================
/// \file SvtxEvaluator.h
/// \brief Compares reconstructed tracks to truth particles
/// \author Michael P. McCumber (revised SVTX version)
//===============================================

#include <fun4all/SubsysReco.h>

#include <string>
#include <TMatrixFfwd.h>
#include <TMatrixT.h>   
#include <TMatrixTUtils.h>

class PHCompositeNode;
class PHTimer;
class TrkrCluster;
class SvtxEvalStack;
class TFile;
class TNtuple;
class SvtxTrack;
class SvtxVertexMap;

//class TrkrClusterContainer;

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
                unsigned int nlayers_tpc = 48,
                unsigned int nlayers_mms = 2);
  ~SvtxEvaluator() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  //  void do_primaries(bool b);

  void set_strict(bool b) { _strict = b; }
  void set_use_initial_vertex(bool use_init_vtx) {_use_initial_vertex = use_init_vtx;}
  void set_use_genfit_vertex(bool use_genfit_vtx) {_use_genfit_vertex = use_genfit_vtx;}
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
  void do_vtx_eval_light(bool b) { _do_vtx_eval_light = b;}
  void scan_for_embedded(bool b) { _scan_for_embedded = b; }
  void scan_for_primaries(bool b) { _scan_for_primaries = b; }
  void set_cluster_version(int value) { m_cluster_version = value; }

 private:
  unsigned int _ievent;
  unsigned int _iseed;
  float m_fSeed;
  // eval stack
  SvtxEvalStack *_svtxevalstack;

  TMatrixF calculateClusterError(TrkrCluster* c, float& clusphi);
  void get_dca(SvtxTrack* track, SvtxVertexMap* vertexmap,
	       float& dca3dxy, float& dca3dz,
	       float& dca3dxysigma, float& dca3dzsigma);
  //TrkrClusterContainer *cluster_map{nullptr};

  //----------------------------------
  // evaluator output ntuples

  bool _strict;
  bool _use_initial_vertex = true;
  bool _use_genfit_vertex = false;
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
  bool _do_vtx_eval_light;
  bool _scan_for_embedded;
  bool _scan_for_primaries;

  unsigned int _nlayers_maps = 3;
  unsigned int _nlayers_intt = 4;
  unsigned int _nlayers_tpc = 48;
  unsigned int _nlayers_mms = 2;

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

  // output subroutines
  void fillOutputNtuples(PHCompositeNode *topNode);  ///< dump the evaluator information into ntuple for external analysis
  void printInputInfo(PHCompositeNode *topNode);     ///< print out the input object information (debugging upstream components)
  void printOutputInfo(PHCompositeNode *topNode);    ///< print out the ancestry information for detailed diagnosis
    int m_cluster_version = 4;
};

#endif  // G4EVAL_SVTXEVALUATOR_H
