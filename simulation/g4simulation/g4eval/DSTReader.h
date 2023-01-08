#ifndef G4EVAL_DSTREADER_H
#define G4EVAL_DSTREADER_H

/*!
 * \file DSTReader.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 * \modified Christof Roland <cer@mit.edu>
 */

// #include "DSTContainerv1.h"

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>


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
class DSTContainerv1;
class TrkrCluster;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;

class DSTReader : public SubsysReco
{
  public:

  //! constructor
  DSTReader( const std::string& = "DSTReader");

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

  // get geant hits associated to a cluster
  using G4HitSet = std::set<PHG4Hit*>;
  G4HitSet find_g4hits( TrkrDefs::cluskey ) const;

  //! get G4Particle id of max contributor to a given track
  std::pair<int,int> get_max_contributor( SvtxTrack* ) const;

  //! get embedded id for given g4track
  int get_embed(PHG4Particle*) const;

  // SvtxTrack recover_track(DSTContainerv1::TrackStruct);

  // TrkrCluster recover_cluster(DSTContainerv1::ClusterStruct);

  //! evaluation node
  DSTContainerv1* m_container = nullptr;

  //! flags
  int m_flags = WriteEvent | WriteClusters | WriteTracks;

  //! hits
  TrkrHitSetContainer* m_hitsetcontainer = nullptr;

  //! clusters
  TrkrClusterContainer* m_cluster_map = nullptr;

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

  // debugging helpers
  void showMe() const;
  void showAll() const;
  void printCluster(TrkrCluster&) const;
};

#endif  // G4EVAL_DSTREADER_H
