// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHSILICONTRACKMATCHING_H
#define PHSILICONTRACKMATCHING_H

#include <fun4all/SubsysReco.h>
#include <trackbase/ActsGeometry.h>

#include <map>
#include <string>

class PHCompositeNode;
class TrackSeedContainer;
class TrackSeed;
class TrkrClusterContainer;

class PHSiliconTrackMatching : public SubsysReco
{
 public:
  PHSiliconTrackMatching(const std::string &name = "PHSiliconTrackMatching");

  ~PHSiliconTrackMatching() override;

  void set_cluster_map_name(const std::string &name)
  {
    _cluster_map_name = name;
  }
  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *) override;

  int End(PHCompositeNode *) override;

  void fieldMap(std::string &fieldmap) { m_fieldMap = fieldmap; }

  void set_silicon_track_map_name(const std::string &map_name) { _silicon_track_map_name = map_name; }
  void set_track_map_name(const std::string &map_name) { _track_map_name = map_name; }
  void SetIteration(int iter) { _n_iteration = iter; }

 private:
  int GetNodes(PHCompositeNode *topNode);

  //  bool _use_old_matching = false;  // normally false

  bool _zero_field = false;     // fit straight lines if true

  TrackSeedContainer *_svtx_seed_map{nullptr};
  TrackSeedContainer *_track_map{nullptr};
  TrackSeedContainer *_track_map_silicon{nullptr};
  TrackSeed *_tracklet_si{nullptr};
  TrkrClusterContainer *_cluster_map{nullptr};
  ActsGeometry *_tGeometry{nullptr};
  int m_event = 0;

  int _n_iteration = 0;
  std::string _track_map_name = "TpcTrackSeedContainer";
  std::string _silicon_track_map_name = "SiliconTrackSeedContainer";
  std::string _cluster_map_name = "TRKR_CLUSTER";
  std::string m_fieldMap = "1.4";
};

#endif  //  PHSiliconTrackMatching_H
