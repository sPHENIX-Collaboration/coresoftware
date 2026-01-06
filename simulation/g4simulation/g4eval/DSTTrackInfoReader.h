#ifndef G4EVAL_DSTTRACKINFOREADER_H
#define G4EVAL_DSTTRACKINFOREADER_H

/*!
 * \file DSTTrackInfoReader.h
 * \author Alex Patton <aopatton@mit.edu>
 */

#include <fun4all/SubsysReco.h>

class PHCompositeNode;
class TrackInfoContainer;

class DSTTrackInfoReader : public SubsysReco
{
 public:
  //! constructor
  DSTTrackInfoReader(const std::string& = "DSTTrackInfoReader");

  //! run initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

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
  int load_nodes(PHCompositeNode *topNode);

  void evaluate_track_info();

  // SvtxTrack recover_track(DSTContainerv3::TrackStruct);

  // TrkrCluster recover_cluster(DSTContainerv3::ClusterStruct);

  TrackInfoContainer *m_track_info_container {nullptr};
  // DSTContainer* m_container = nullptr;

  //! flags
  int m_flags = WriteEvent | WriteClusters | WriteTracks;

  //@}

  // debugging helpers
  // bool dryrun = false;
  // bool generateKey = false;
};

#endif  // G4EVAL_DSTTRACKINFOREADER_H
