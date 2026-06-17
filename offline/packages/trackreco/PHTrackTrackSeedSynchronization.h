#ifndef PHTRACKTRACKSEEDSYNCHRONIZATION_H
#define PHTRACKTRACKSEEDSYNCHRONIZATION_H

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

class PHTrackTrackSeedSynchronization : public SubsysReco
{
 public:

  //! constructor
  PHTrackTrackSeedSynchronization(const std::string& = "PHTrackTrackSeedSynchronization");

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *) override;

  int End(PHCompositeNode *) override;

  //! input silicon seeds map name
  void set_si_seed_map_name(const std::string &map_name) { _si_seed_map_name = map_name; }

  //! input tpc seeds map name
  void set_tpc_seed_map_name(const std::string &map_name) { _tpc_seed_map_name = map_name; }

  //! input track map name
  void set_svtx_track_map_name(const std::string &map_name) { _svtx_track_map_name = map_name; }

  private:

  int GetNodes(PHCompositeNode*);

  /// make sure that the track seed pointers stored in track correspond to those store in the seed containers
  bool synchronize_track(SvtxTrack*) const;

  TrackSeedContainer *_tpc_seed_map{nullptr};
  TrackSeedContainer *_si_seed_map{nullptr};
  SvtxTrackMap *_svtx_track_map{nullptr};

  std::string _tpc_seed_map_name = "TpcTrackSeedContainer";
  std::string _si_seed_map_name = "SiliconTrackSeedContainer";
  std::string _svtx_track_map_name = "SvtxTrackMap";

};

#endif  //  PHTrackTrackSeedSynchronization_H
