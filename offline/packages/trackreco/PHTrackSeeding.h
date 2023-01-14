/*!
 *  \file		  PHTrackSeeding.h
 *  \brief		Base class for track seeding
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef TRACKRECO_PHTRACKSEEDING_H
#define TRACKRECO_PHTRACKSEEDING_H

// PHENIX includes
#include <fun4all/SubsysReco.h>

// STL includes
#include <string>

// forward declarations
class PHCompositeNode;

//class SvtxClusterMap;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitSetContainer;
class TrkrClusterIterationMapv1;
class SvtxVertexMap;
class TrackSeedContainer;


/// \class PHTrackSeeding
///
/// \brief Base class for track seeding
///
class PHTrackSeeding : public SubsysReco
{
 public:
  PHTrackSeeding(const std::string &name = "PHTrackSeeding");
  ~PHTrackSeeding() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void set_track_map_name(const std::string &map_name) { _track_map_name = map_name; }
  void set_do_hit_association(bool do_assoc){do_hit_assoc = do_assoc;}
  void SetUseTruthClusters(bool setit){_use_truth_clusters = setit;}
  void SetIteration(int iter){_n_iteration = iter;}
 protected:
  /// setup interface for trackers, called in InitRun, setup things like pointers to nodes.
  /// overrided in derived classes
  virtual int Setup(PHCompositeNode *topNode);

  /// process event interface for trackers, called in process_event.
  /// implemented in derived classes
  virtual int Process(PHCompositeNode *topNode) = 0;

  /// Called in SubsysReco::End
  virtual int End() = 0;

  TrkrClusterContainer *_cluster_map = nullptr;
  TrkrClusterHitAssoc *_cluster_hit_map = nullptr;
  TrkrClusterIterationMapv1* _iteration_map;
  int _n_iteration;
  bool do_hit_assoc = true;
  SvtxVertexMap *_vertex_map = nullptr;
  TrackSeedContainer *_track_map = nullptr;
  TrkrHitSetContainer  *_hitsets = nullptr;

  std::string _track_map_name = "TpcTrackSeedContainer";

  bool _use_truth_clusters = false;

 private:

  /// create new node output pointers
  int CreateNodes(PHCompositeNode *topNode);

  /// fetch node pointers
  int GetNodes(PHCompositeNode *topNode);
};

#endif
