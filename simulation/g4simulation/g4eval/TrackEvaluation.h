#ifndef G4EVAL_TRACKEVALUATION_H
#define G4EVAL_TRACKEVALUATION_H

/*!
 * \file TrackEvaluation.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TrackEvaluationContainerv1.h"

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/ClusterErrorPara.h>

#include <trackbase_historic/SvtxTrackState.h>
#include <map>
#include <set>
#include <string>
#include <vector>

class ActsGeometry;
class PHG4TpcCylinderGeomContainer;
class PHG4CylinderGeomContainer;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class PHG4TruthInfoContainer;
class SvtxTrack;
class SvtxTrackMap;
class TrkrCluster;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;

class TrackEvaluation : public SubsysReco
{
  public:

  //! constructor
  TrackEvaluation( const std::string& = "TrackEvaluation" );

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
    EvalEvent = 1<<0,
    EvalClusters = 1<<1,
    EvalTracks = 1<<2
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
  G4HitSet find_g4hits( TrkrDefs::cluskey, int id ) const;

  //! get G4Particle id of max contributor to a given track
  std::pair<int,int> get_max_contributor( SvtxTrack* ) const;

  //! get embedded id for given g4track
  int get_embed(PHG4Particle*) const;

  //! create cluster structure from cluster
  TrackEvaluationContainerv1::ClusterStruct create_cluster( TrkrDefs::cluskey, TrkrCluster*,SvtxTrack* ) const;

  //! add track information to a cluster
  void add_trk_information( TrackEvaluationContainerv1::ClusterStruct&, SvtxTrackState* ) const;

  //! add track information to a cluster for the micromegas case
  /*!
   * the difference between this and the generic method is that the track state to
   * the tiles detector plane, and not to the same radius as the cluster
   */
  void add_trk_information_micromegas( TrackEvaluationContainerv1::ClusterStruct&, int /* tileid */, SvtxTrackState* ) const;

  // add truth information
  void add_truth_information( TrackEvaluationContainerv1::ClusterStruct&, std::set<PHG4Hit*> ) const;

  // add truth information
  /*!
   * the difference between this and the generic method is that the track state to
   * the tiles detector plane, and not to the same radius as the cluster
   */
  void add_truth_information_micromegas( TrackEvaluationContainerv1::ClusterStruct&, int /* tileid */, std::set<PHG4Hit*> ) const;

  //! evaluation node
  TrackEvaluationContainerv1* m_container = nullptr;

  //! flags
  int m_flags = EvalEvent | EvalClusters | EvalTracks;

  /// Acts tracking geometry for surface lookup
  ActsGeometry *m_tGeometry = nullptr;

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

  //! tpc geometry
  PHG4TpcCylinderGeomContainer* m_tpc_geom_container = nullptr;

  //! micromegas geometry
  PHG4CylinderGeomContainer* m_micromegas_geom_container = nullptr;

  // map cluster keys to g4hits
  using G4HitMap = std::map<TrkrDefs::cluskey,G4HitSet>;
  mutable G4HitMap m_g4hit_map;
  ClusterErrorPara _ClusErrPara;

};

#endif  // G4EVAL_TrackEvaluation_H
