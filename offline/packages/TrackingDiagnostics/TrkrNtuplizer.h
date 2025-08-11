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
#include <queue>

class CDBInterface;
class CDBTTree;
class PHCompositeNode;
class PHTimer;
class TrkrCluster;
class TFile;
class TNtuple;
class SvtxTrack;
class TrackSeed;
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
  void do_dedx_calib(bool b) { _do_dedx_calib = b; }
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
  static std::vector<TrkrDefs::cluskey> get_track_ckeys(SvtxTrack *track);
  void segment(const int seg) { m_segment = seg; }
  void runnumber(const int run) { m_runnumber = run; }
  void job(const int job) { m_job = job; }

 private:
  struct fee_info
  {
    unsigned int fee = std::numeric_limits<unsigned int>::max();
    unsigned int channel = std::numeric_limits<unsigned int>::max();
    unsigned int sampa = std::numeric_limits<unsigned int>::max();
  };
  int m_segment = 0;
  int m_runnumber = 0;
  int m_job = 0;

  unsigned int _ievent{0};
  unsigned int _iseed{0};
  float m_fSeed{std::numeric_limits<float>::quiet_NaN()};
  // eval stack

  float calc_dedx(TrackSeed *tpcseed);
  TF1 *f_pion_plus{nullptr};
  TF1 *f_kaon_plus{nullptr};
  TF1 *f_proton_plus{nullptr};
  TF1 *f_pion_minus{nullptr};
  TF1 *f_kaon_minus{nullptr};
  TF1 *f_proton_minus{nullptr};
  float  dedxcorr[2][12][3]{};
  float get_n1pix(TrackSeed *tpcseed);

  static TMatrixF calculateClusterError(TrkrCluster *c, float &clusphi);
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
  bool _do_dedx_calib{false};
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
  float m_ZDC_coincidence{0};
  float m_mbd_rate{0};
  float m_rawzdc{0};
  float m_livezdc{0};
  float m_scaledzdc{0};
  float m_rawmbd{0};
  float m_livembd{0};
  float m_scaledmbd{0};
  float m_rawmbdv10{0};
  float m_livembdv10{0};
  float m_scaledmbdv10{0};
  float m_mbd_rate1{0};
  float m_rawzdc1{0};
  float m_livezdc1{0};
  float m_scaledzdc1{0};
  float m_rawmbd1{0};
  float m_livembd1{0};
  float m_scaledmbd1{0};
  float m_rawmbdv101{0};
  float m_livembdv101{0};
  float m_scaledmbdv101{0};
  float m_rawzdclast{0};
  float m_rawmbdlast{0};
  float m_rawmbdv10last{0};
  uint64_t m_bcolast = std::numeric_limits<uint64_t>::quiet_NaN();

  std::queue<int> m_rawzdc_hist;
  std::queue<int> m_rawmbd_hist;
  std::queue<int> m_rawmbdv10_hist;
  std::queue<uint64_t> m_bco_hist;

  uint64_t m_bco = std::numeric_limits<uint64_t>::quiet_NaN();
  uint64_t m_bco1 = std::numeric_limits<uint64_t>::quiet_NaN();
  uint64_t m_bcotr = std::numeric_limits<uint64_t>::quiet_NaN();
  uint64_t m_bcotr1 = std::numeric_limits<uint64_t>::quiet_NaN();
  float m_totalmbd = std::numeric_limits<float>::quiet_NaN();
  std::vector<int> m_firedTriggers;
  uint64_t m_gl1BunchCrossing = std::numeric_limits<uint64_t>::quiet_NaN();

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

  CDBTTree *m_cdbttree{nullptr};
  CDBInterface *m_cdb{nullptr};
  
  int mc_sectors[12]{5, 4, 3, 2, 1, 0, 11, 10, 9, 8, 7, 6};
  int FEE_map[26]{4, 5, 0, 2, 1, 11, 9, 10, 8, 7, 6, 0, 1, 3, 7, 6, 5, 4, 3, 2, 0, 2, 1, 3, 5, 4};
  int FEE_R[26]{2, 2, 1, 1, 1, 3, 3, 3, 3, 3, 3, 2, 2, 1, 2, 2, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3};

  std::map<TrkrDefs::cluskey, fee_info> fee_map;
};

#endif  // G4EVAL_SVTXEVALUATOR_H
