#ifndef G4EVAL_DSTTRACKINFOWRITER_H
#define G4EVAL_DSTTRACKINFOWRITER_H

/*!
 * \file DSTTrackInfoWriter.h
 * \author Alex Patton <aopatton@mit.edu>
 *
 */

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrackInfo_v1.h>
#include <trackbase_historic/SvtxTrack_v4.h>
#include <trackbase_historic/TrackInfoContainer_v1.h>
#include <trackbase_historic/TrackStateInfo_v1.h>

#include <map>
#include <set>
#include <string>
#include <vector>

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>

class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class PHG4TruthInfoContainer;
class SvtxTrack;
class SvtxTrackMap;
// class DSTContainerv3;
// class DSTContainer;
class TrkrCluster;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;
class TrackSeedContainer;

class DSTTrackInfoWriter : public SubsysReco
{
 public:
  //! constructor
  DSTTrackInfoWriter(const std::string& = "DSTTrackInfoWriter");

  //! global initialization
  int InitRun(PHCompositeNode* topNode) override;

  //! event processing
  int process_event(PHCompositeNode* topNode) override;

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

  //! evaluate track info
  void evaluate_track_info();

  TrackInfoContainer_v1* m_track_info_container = nullptr;

  //! flags
  int m_flags = WriteEvent | WriteClusters | WriteTracks;

  //! tracks
  SvtxTrackMap* m_track_map = nullptr;
};

#endif  // G4EVAL_DSTTrackInfoWRITER_H
