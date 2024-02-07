#ifndef G4EVAL_DSTTRACKINFOREADER_H
#define G4EVAL_DSTTRACKINFOREADER_H

/*!
 * \file DSTTrackInfoReader.h
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
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;

class DSTTrackInfoReader : public SubsysReco
{
 public:
  //! constructor
  DSTTrackInfoReader(const std::string& = "DSTTrackInfoReader");

  //! run initialization
  int InitRun(PHCompositeNode*) override;

  //! event processing
  int process_event(PHCompositeNode*) override;

  enum Flags
  {
    WriteEvent = 1 << 0,
    WriteClusters = 1 << 1,
    WriteTracks = 1 << 2
  };

  //! set flags. Should be a bitwise or of Flags enum
  void set_flags(int flags)
  {
    m_flags = flags;
  }

 private:
  //! load nodes
  int load_nodes(PHCompositeNode*);

  void evaluate_track_info();

  // SvtxTrack recover_track(DSTContainerv3::TrackStruct);

  // TrkrCluster recover_cluster(DSTContainerv3::ClusterStruct);

  TrackInfoContainer_v1* m_track_info_container = nullptr;
  // DSTContainer* m_container = nullptr;

  //! flags
  int m_flags = WriteEvent | WriteClusters | WriteTracks;

  //@}

  // debugging helpers
  // bool dryrun = false;
  // bool generateKey = false;
};

#endif  // G4EVAL_DSTTRACKINFOREADER_H
