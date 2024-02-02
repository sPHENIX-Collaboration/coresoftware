#ifndef G4EVAL_DSTWRITER_H
#define G4EVAL_DSTWRITER_H

/*!
 * \file DSTWriter.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 * \modified Christof Roland <cer@mit.edu>
 */

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterContainerv5.h>
#include <trackbase_historic/SvtxTrack_v4.h>
#include <trackbase_historic/SvtxTrackArray_v1.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/SvtxTrackInfo_v1.h>
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
class DSTContainerv3;
//class DSTContainer;
class TrkrCluster;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;
class TrackSeedContainer;

class DSTWriter : public SubsysReco
{
  public:

  //! constructor
  DSTWriter( const std::string& = "DSTWriter");

  //! global initialization
  int Init(PHCompositeNode*) override;

  //! run initialization
  int InitRun(PHCompositeNode*) override;

  //! event processing
  int process_event(PHCompositeNode*) override;

  //! end of processing
  int End(PHCompositeNode*) override;

  enum Flags
  {
    WriteEvent = 1<<0,
    WriteClusters = 1<<1,
    WriteTracks = 1<<2
  };

  //! set flags. Should be a bitwise or of Flags enum
  void set_flags( int flags )
  { m_flags = flags; }

  private:

  //! load nodes
  int load_nodes( PHCompositeNode* );

  //! evaluate event
  void evaluate_event();

  //! evaluate clusters
  void evaluate_clusters();

  //! evaluate tracks
  void evaluate_tracks();

  //! evaluate both
  void evaluate_track_and_cluster();

  //! evaluate track info
  void evaluate_track_info();

  // get geant hits associated to a cluster
  using G4HitSet = std::set<PHG4Hit*>;
  G4HitSet find_g4hits( TrkrDefs::cluskey ) const;

  //! get G4Particle id of max contributor to a given track
  std::pair<int,int> get_max_contributor( SvtxTrack* ) const;

  //! get embedded id for given g4track
  int get_embed(PHG4Particle*) const;

  //! evaluation node
  DSTContainerv3* m_container = nullptr;
  //DSTContainer* m_container = nullptr;

  TrackInfoContainer_v1* m_track_info_container = nullptr;

  //! flags
  int m_flags = WriteEvent | WriteClusters | WriteTracks;

  TrackSeedContainer* m_svtxtrackseed = nullptr;
  TrackSeedContainer *_tpc_seeds = nullptr;

  //! hits
  TrkrHitSetContainer* m_hitsetcontainer = nullptr;

  //! clusters
  TrkrClusterContainer* m_cluster_map = nullptr;

  //! clusters in array
  TrkrClusterContainerv5* m_cluster_map_arr = nullptr;

  //! cluster to hit association
  TrkrClusterHitAssoc* m_cluster_hit_map = nullptr;

  //! hit to truth association
  TrkrHitTruthAssoc* m_hit_truth_map = nullptr;

  //! tracks
  SvtxTrackMap* m_track_map = nullptr;

  //!@name geant4 hits
  //@{
  PHG4HitContainer* m_g4hits_tpc = nullptr;
  PHG4HitContainer* m_g4hits_intt = nullptr;
  PHG4HitContainer* m_g4hits_mvtx = nullptr;
  PHG4HitContainer* m_g4hits_micromegas = nullptr;
  //@}

  //! truth information
  PHG4TruthInfoContainer* m_g4truthinfo = nullptr;

  // map cluster keys to g4hits
  using G4HitMap = std::map<TrkrDefs::cluskey,G4HitSet>;
  mutable G4HitMap m_g4hit_map;

  TFile* fcl = nullptr;
  TTree* tcl = nullptr;
  TClonesArray* arrEvt = nullptr;
  //TClonesArray* arrTrk = nullptr;
  TClonesArray* arrCls = nullptr;

};

#endif  // G4EVAL_DSTWRITER_H
