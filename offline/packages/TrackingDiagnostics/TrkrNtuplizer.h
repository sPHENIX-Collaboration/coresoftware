#ifndef TRACKRECO_TRKRNTUPLIZER_H
#define TRACKRECO_TRKRNTUPLIZER_H

//===============================================
/// \file SvtxEvaluator.h
/// \brief Compares reconstructed tracks to truth particles
/// \author Michael P. McCumber (revised SVTX version)
//===============================================

#include <fun4all/SubsysReco.h>

#include <TMatrixFfwd.h>
#include <TMatrixT.h>
#include <TMatrixTUtils.h>
#include <string>

class PHCompositeNode;
class PHTimer;
class TrkrCluster;
class TFile;
class TNtuple;
class SvtxTrack;
class SvtxTrackEval;
class SvtxVertexMap;

// class TrkrClusterContainer;

/// \class SvtxEvaluator
///
/// \brief Compares reconstructed tracks to truth particles
///
/// Plan: This module will trace the reconstructed tracks back to
/// the greatest contributor Monte Carlo particle and then
/// test one against the other.
///
class TrkrNtuplizer : public SubsysReco
{
 public:
  TrkrNtuplizer(const std::string &name = "TRKRNTUPLEIZER",
                const std::string &filename = "trkrntuple.root",
                const std::string &trackmapname = "SvtxTrackMap",
                unsigned int nlayers_maps = 3,
                unsigned int nlayers_intt = 8,
                unsigned int nlayers_tpc = 48,
                unsigned int nlayers_mms = 2);
  ~TrkrNtuplizer() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  //  void do_primaries(bool b);

  void do_info_eval(bool b) { _do_info_eval = b; }
  void do_vertex_eval(bool b) { _do_vertex_eval = b; }
  void do_hit_eval(bool b) { _do_hit_eval = b; }
  void do_cluster_eval(bool b) { _do_cluster_eval = b; }
  void do_track_eval(bool b) { _do_track_eval = b; }
  void do_tpcseed_eval(bool b) { _do_tpcseed_eval = b; }
  void do_siseed_eval(bool b) { _do_siseed_eval = b; }

  void set_cluster_version(int value) { m_cluster_version = value; }

 private:
  unsigned int _ievent;
  unsigned int _iseed;
  float m_fSeed;
  // eval stack

  TMatrixF calculateClusterError(TrkrCluster *c, float &clusphi);
  void get_dca(SvtxTrack *track, SvtxVertexMap *vertexmap,
               float &dca3dxy, float &dca3dz,
               float &dca3dxysigma, float &dca3dzsigma);
  // TrkrClusterContainer *cluster_map{nullptr};

  //----------------------------------
  // evaluator output ntuples

  bool _do_info_eval;
  bool _do_vertex_eval;
  bool _do_hit_eval;
  bool _do_cluster_eval;
  bool _do_track_eval;
  bool _do_tpcseed_eval;
  bool _do_siseed_eval;

  unsigned int _nlayers_maps = 3;
  unsigned int _nlayers_intt = 4;
  unsigned int _nlayers_tpc = 48;
  unsigned int _nlayers_mms = 2;

  TNtuple *_ntp_info;
  TNtuple *_ntp_vertex;
  TNtuple *_ntp_hit;
  TNtuple *_ntp_cluster;
  TNtuple *_ntp_track;
  TNtuple *_ntp_tpcseed;
  TNtuple *_ntp_siseed;

  // evaluator output file
  std::string _filename;
  // Track map name
  std::string _trackmapname;
  TFile *_tfile;

  PHTimer *_timer;

  // output subroutines
  void fillOutputNtuples(PHCompositeNode *topNode);  ///< dump the evaluator information into ntuple for external analysis
  void printInputInfo(PHCompositeNode *topNode);     ///< print out the input object information (debugging upstream components)
  void printOutputInfo(PHCompositeNode *topNode);    ///< print out the ancestry information for detailed diagnosis
  double AdcClockPeriod = 53.0;                      // ns
  int m_cluster_version = 4;
  SvtxTrackEval *_TrackEval;
};

#endif  // G4EVAL_SVTXEVALUATOR_H
