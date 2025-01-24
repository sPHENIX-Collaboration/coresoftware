#ifndef DSTCLUSTERPRUNING_H
#define DSTCLUSTERPRUNING_H

/*!
 * \file DSTClusterPruning.h
 * \author Alex Patton <aopatton@mit.edu>
 */

// #include "DSTContainerv3.h"

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterv5.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/TrackInfoContainer_v1.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v1.h>

#include <map>
#include <set>
#include <string>
#include <vector>

class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class PHG4TruthInfoContainer;
class SvtxTrack;
class SvtxTrackMap;
// class DSTContainer;
class TrkrCluster;
class TrkrClusterv5;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;
class TrackSeedContainer;

class DSTClusterPruning : public SubsysReco
{
 public:
  //! constructor
  DSTClusterPruning(const std::string& = "DSTClusterPruning");

  //! run initialization
  int Init(PHCompositeNode*) override;
  //int InitRun(PHCompositeNode*) override;

  //! event processing
  int process_event(PHCompositeNode*) override;

  //! end of processing
  //int End(PHCompositeNode*) override;

 private:
  //! load nodes
  int load_nodes(PHCompositeNode*);

  void prune_clusters();
  void fill_clusters();
  void print_clusters();
  SvtxTrackMap* m_track_map = nullptr;

  TrkrClusterv5* m_cluster = nullptr;
  TrkrClusterContainer* m_cluster_map = nullptr;
  TrkrClusterContainer* m_reduced_cluster_map = nullptr;
  TrackSeedContainer* m_track_seed_container = nullptr;
  TrackSeedContainer* m_tpc_track_seed_container = nullptr;
  TrackSeedContainer* m_silicon_track_seed_container = nullptr;

  //@}

  // debugging helpers
  // bool dryrun = false;
  // bool generateKey = false;
};

#endif  // DSTCLUSTERPRUNING_H