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
  PHTrackPruner(const std::string &name = "PHTrackPruner");

  ~PHTrackPruner() override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *) override;

  int End(PHCompositeNode *) override;

  void set_pruned_svtx_seed_map_name(const std::string &map_name) { _pruned_svtx_seed_map_name = map_name; }
  void set_svtx_seed_map_name(const std::string &map_name) { _svtx_seed_map_name = map_name; }
  void set_si_seed_map_name(const std::string &map_name) { _si_seed_map_name = map_name; }
  void set_tpc_seed_map_name(const std::string &map_name) { _tpc_seed_map_name = map_name; }
  void set_svtx_track_map_name(const std::string &map_name) { _svtx_track_map_name = map_name; }

  void set_track_pt_low_cut(const double val) { m_track_pt_low_cut = val; }
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
  TrackSeedContainer *_svtx_seed_map{nullptr};
  TrackSeedContainer *_tpc_seed_map{nullptr};
  TrackSeedContainer *_si_seed_map{nullptr};
  TrackSeed *_tpc_seed{nullptr};
  TrackSeed *_si_seed{nullptr};
  SvtxTrackMap *_svtx_track_map{nullptr};
  SvtxTrack *_svtx_track{nullptr};
  TrkrClusterContainer *_cluster_map{nullptr};
  ActsGeometry *_tGeometry{nullptr};
  int m_event = 0;

  std::string _tpc_seed_map_name = "TpcTrackSeedContainer";
  std::string _si_seed_map_name = "SiliconTrackSeedContainer";
  std::string _svtx_seed_map_name = "SvtxTrackSeedContainer";
  std::string _pruned_svtx_seed_map_name = "PrunedSvtxTrackSeedContainer";
  std::string _svtx_track_map_name = "SvtxTrackMap";

  double m_track_pt_low_cut = 0.5;
  double m_track_quality_high_cut = 100;
  int m_nmvtx_clus_low_cut = 3;
  int m_nintt_clus_low_cut = 2;
  int m_ntpc_clus_low_cut = 35;
  int m_ntpot_clus_low_cut = 2;

  int m_nmvtx_states_low_cut = 3;
  int m_nintt_states_low_cut = 2;
  int m_ntpc_states_low_cut = 35;
  int m_ntpot_states_low_cut = 2;

};

#endif  //  PHTRACKPRUNER_H
