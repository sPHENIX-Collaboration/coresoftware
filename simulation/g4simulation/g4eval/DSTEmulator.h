#ifndef G4EVAL_DSTEMULATOR_H
#define G4EVAL_DSTEMULATOR_H

/*!
 * \file TrackEvaluation.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase_historic/ActsTransformations.h>

#include <map>
#include <set>
#include <string>
#include <vector>
#include <TRandom.h>

#include "DSTCompressor.h"

class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class PHG4TruthInfoContainer;
class SvtxTrack;
class SvtxTrackMap;
class TrackEvaluationContainerv1;
class TrkrCluster;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;

class TFile;
class TNtuple;

class DSTEmulator : public SubsysReco
{
  public:

  //! constructor
  DSTEmulator( const std::string& = "DSTEmulator",
               const std::string &filename = "DSTana.root",
               int nBits = 8,
               int sabotage = 0,
               bool compress = true);

  //! global initialization
  int Init(PHCompositeNode*) override;

  //! run initialization
  int InitRun(PHCompositeNode*) override;

  //! event processing
  int process_event(PHCompositeNode*) override;

  //! end of processing
  int End(PHCompositeNode*) override;

  private:

  //! load nodes
  int load_nodes( PHCompositeNode* );

  //! evaluate tracks
  void evaluate_tracks();
  float compress_dx(float in_val);
  float compress_dy(float in_val);
  // get geant hits associated to a cluster
  using G4HitSet = std::set<PHG4Hit*>;
  G4HitSet find_g4hits( TrkrDefs::cluskey ) const;
  Acts::Vector3 getGlobalPosition(TrkrDefs::cluskey, TrkrCluster*) const;
  //! get G4Particle id of max contributor to a given track
  std::pair<int,int> get_max_contributor( SvtxTrack* ) const;

  //! get embedded id for given g4track
  int get_embed(PHG4Particle*) const;

  //! evaluation node
  TrackEvaluationContainerv1* m_container = nullptr;

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

  ActsGeometry *m_tGeometry = nullptr;
  ActsTransformations m_transform;

  TNtuple *_dst_data;

  // output file
  std::string _filename;
  TFile *_tfile;

  DSTCompressor* m_compressor;

  // Number of bits for the integer representation after compression
  int nBits = 8;
  // replace the decompressed residuals by a large number
  int sabotage = 0;
  // random seed
  TRandom rnd;
  // switch to apply the compressed residuals to cluster positions
  bool apply_compression = true;

};

#endif  // G4EVAL_DSTEMULATOR_H
