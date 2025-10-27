#ifndef G4EVAL_DSTTRACKINFOWRITER_H
#define G4EVAL_DSTTRACKINFOWRITER_H

/*!
 * \file DSTTrackInfoWriter.h
 * \author Alex Patton <aopatton@mit.edu>
 *
 */

#include <fun4all/SubsysReco.h>

class PHCompositeNode;
class SvtxTrackMap;
class TrackInfoContainer;

class DSTTrackInfoWriter : public SubsysReco
{
 public:
  //! constructor
  DSTTrackInfoWriter(const std::string& = "DSTTrackInfoWriter");

  //! global initialization
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

  //! evaluate track info
  void evaluate_track_info();

  TrackInfoContainer *m_track_info_container {nullptr};

  //! flags
  int m_flags {WriteEvent | WriteClusters | WriteTracks};

  //! tracks
  SvtxTrackMap *m_track_map {nullptr};
};

#endif  // G4EVAL_DSTTrackInfoWRITER_H
