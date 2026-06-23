// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHTRACKPRUNER_H
#define PHTRACKPRUNER_H

#include <fun4all/SubsysReco.h>
#include <phparameter/PHParameterInterface.h>
#include <trackbase/ActsGeometry.h>

#include <map>
#include <string>

class PHCompositeNode;
class TrackSeedContainer;
class TrackSeed;
class SvtxTrackSeed;
class SvtxTrackMap;
class SvtxTrack;
class TrkrClusterContainer;
class TF1;
class TFile;
class TNtuple;

class PHTrackPruner : public SubsysReco
{
 public:


  //! constructor
  PHTrackPruner(const std::string &name = "PHTrackPruner");

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *) override;

  int End(PHCompositeNode *) override;

  //! input cluster map name. Default is TRKR_CLUSTER
  void set_cluster_map_name(const std::string &map_name) { _cluster_map_name = map_name; }

  //! input silicon seeds map name
  void set_si_seed_map_name(const std::string &map_name) { _si_seed_map_name = map_name; }

  //! input tpc seeds map name
  void set_tpc_seed_map_name(const std::string &map_name) { _tpc_seed_map_name = map_name; }

  //! input track map name
  void set_svtx_track_map_name(const std::string &map_name) { _svtx_track_map_name = map_name; }

  //! output pruned track map name
  void set_pruned_svtx_seed_map_name(const std::string &map_name) { _pruned_svtx_seed_map_name = map_name; }

  /// low pt cut
  void set_track_pt_low_cut(const double val) { m_track_pt_low_cut = val; }

  /// high pt cut.
  /** enforced only if >0 */
  void set_track_pt_high_cut(const double val) { m_track_pt_high_cut = val; }

  void set_track_quality_high_cut(const double val) { m_track_quality_high_cut = val; }

  void set_nmvtx_clus_low_cut(const int n) { m_nmvtx_clus_low_cut = n; }
  void set_nintt_clus_low_cut(const int n) { m_nintt_clus_low_cut = n; }
  void set_ntpc_clus_low_cut(const int n) { m_ntpc_clus_low_cut = n; }
  void set_ntpot_clus_low_cut(const int n) { m_ntpot_clus_low_cut = n; }

  void set_nmvtx_states_low_cut(const int n) { m_nmvtx_states_low_cut = n; }
  void set_nintt_states_low_cut(const int n) { m_nintt_states_low_cut = n; }
  void set_ntpc_states_low_cut(const int n) { m_ntpc_states_low_cut = n; }
  void set_ntpot_states_low_cut(const int n) { m_ntpot_states_low_cut = n; }

 private:
  int GetNodes(PHCompositeNode *topNode);
  bool checkTrack(SvtxTrack *track);

  short int findCrossingGeometrically(unsigned int tpcid, unsigned int si_id);
  double getBunchCrossing(unsigned int trid, double z_mismatch);

  TrackSeedContainer *_pruned_svtx_seed_map{nullptr};
  TrackSeedContainer *_tpc_seed_map{nullptr};
  TrackSeedContainer *_si_seed_map{nullptr};

  SvtxTrackMap *_svtx_track_map{nullptr};
  TrkrClusterContainer *_cluster_map{nullptr};
  ActsGeometry *_tGeometry{nullptr};

  std::string _cluster_map_name = "TRKR_CLUSTER";
  std::string _tpc_seed_map_name = "TpcTrackSeedContainer";
  std::string _si_seed_map_name = "SiliconTrackSeedContainer";
  std::string _pruned_svtx_seed_map_name = "PrunedSvtxTrackSeedContainer";
  std::string _svtx_track_map_name = "SvtxTrackMap";

  /// low pt cut
  double m_track_pt_low_cut = 0.5;

  /// high pt cut.
  /** enforced only if >0 */
  double m_track_pt_high_cut = 0.;

  double m_track_quality_high_cut = 100;
  int m_nmvtx_clus_low_cut = 3;
  int m_nintt_clus_low_cut = 2;
  int m_ntpc_clus_low_cut = 35;
  int m_ntpot_clus_low_cut = 2;

  int m_nmvtx_states_low_cut = 3;
  int m_nintt_states_low_cut = 2;
  int m_ntpc_states_low_cut = 35;
  int m_ntpot_states_low_cut = 2;

  //! keep track of track/seed statistics
  unsigned long m_total_tracks = 0;
  unsigned long m_accepted_tracks = 0;

};

#endif  //  PHTRACKPRUNER_H
