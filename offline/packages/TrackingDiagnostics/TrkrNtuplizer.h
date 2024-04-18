#ifndef TRACKRECO_TRKRNTUPLIZER_H
#define TRACKRECO_TRKRNTUPLIZER_H

//===============================================
/// \file SvtxEvaluator.h
/// \brief Compares reconstructed tracks to truth particles
/// \author Michael P. McCumber (revised SVTX version)
//===============================================

#include <fun4all/SubsysReco.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrkrDefs.h>

#include <TMatrixFfwd.h>
#include <TMatrixT.h>
#include <TMatrixTUtils.h>

#include <limits>
#include <map>
#include <set>
#include <string>

class PHCompositeNode;
class PHTimer;
class TrkrCluster;
class TFile;
class TNtuple;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class TrkrClusterContainer;
class ActsGeometry;
class PHG4TpcCylinderGeomContainer;
class GlobalVertexMap;

// class ClusterErrorPara;

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
                unsigned int nlayers_intt = 4,
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
  void do_clus_trk_eval(bool b) { _do_clus_trk_eval = b; }
  void do_track_eval(bool b) { _do_track_eval = b; }
  void do_tpcseed_eval(bool b) { _do_tpcseed_eval = b; }
  void do_siseed_eval(bool b) { _do_siseed_eval = b; }
  void set_first_event(int value) { _ievent = value; }
  void set_trkclus_seed_container(const std::string &name)
  {
    _clustrackseedcontainer = name;
  }
  SvtxTrack *best_track_from(TrkrDefs::cluskey cluster_key);
  std::set<SvtxTrack *> all_tracks_from(TrkrDefs::cluskey cluster_key);
  void create_cache_track_from_cluster();
  std::vector<TrkrDefs::cluskey> get_track_ckeys(SvtxTrack *track);

 private:
  unsigned int _ievent{0};
  unsigned int _iseed{0};
  float m_fSeed{std::numeric_limits<float>::quiet_NaN()};
  // eval stack

  TMatrixF calculateClusterError(TrkrCluster *c, float &clusphi);
  void get_dca(SvtxTrack *track, SvtxVertexMap *vertexmap,
               float &dca3dxy, float &dca3dz,
               float &dca3dxysigma, float &dca3dzsigma);
  // TrkrClusterContainer *cluster_map{nullptr};

  void FillCluster(Float_t fXcluster[30], TrkrDefs::cluskey cluster_key);
  void FillTrack(Float_t fXcluster[30], SvtxTrack *track, GlobalVertexMap *vertexmap);
  //----------------------------------
  // evaluator output ntuples

  bool _do_info_eval{true};
  bool _do_vertex_eval{true};
  bool _do_hit_eval{true};
  bool _do_cluster_eval{true};
  bool _do_clus_trk_eval{true};
  bool _do_track_eval{true};
  bool _do_tpcseed_eval{false};
  bool _do_siseed_eval{false};

  unsigned int _nlayers_maps{3};
  unsigned int _nlayers_intt{4};
  unsigned int _nlayers_tpc{48};
  unsigned int _nlayers_mms{2};

  TNtuple *_ntp_info{nullptr};
  TNtuple *_ntp_vertex{nullptr};
  TNtuple *_ntp_hit{nullptr};
  TNtuple *_ntp_cluster{nullptr};
  TNtuple *_ntp_clus_trk{nullptr};
  TNtuple *_ntp_track{nullptr};
  TNtuple *_ntp_tpcseed{nullptr};
  TNtuple *_ntp_siseed{nullptr};

  // evaluator output file
  std::string _filename;
  // Track map name
  std::string _trackmapname;
  ClusterErrorPara _ClusErrPara;
  TrkrClusterContainer *_cluster_map{nullptr};
  SvtxTrackMap *_trackmap{nullptr};
  ActsGeometry *_tgeometry{nullptr};
  PHG4TpcCylinderGeomContainer *_geom_container{nullptr};

  std::string _clustrackseedcontainer = "TpcTrackSeedContainer";

  TFile *_tfile{nullptr};
  PHTimer *_timer{nullptr};

  // output subroutines
  void fillOutputNtuples(PHCompositeNode *topNode);  ///< dump the evaluator information into ntuple for external analysis
  void printInputInfo(PHCompositeNode *topNode);     ///< print out the input object information (debugging upstream components)
  void printOutputInfo(PHCompositeNode *topNode);    ///< print out the ancestry information for detailed diagnosis
  double AdcClockPeriod = 53.0;                      // ns

  bool _cache_track_from_cluster_exists = false;
  std::map<TrkrDefs::cluskey, std::set<SvtxTrack *> > _cache_all_tracks_from_cluster;
  std::map<TrkrDefs::cluskey, SvtxTrack *> _cache_best_track_from_cluster;
};

#endif  // G4EVAL_SVTXEVALUATOR_H
