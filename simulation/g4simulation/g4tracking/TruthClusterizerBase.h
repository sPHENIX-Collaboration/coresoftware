#ifndef G4TRACKING_TRUTHCLUSTERIZERBASE
#define G4TRACKING_TRUTHCLUSTERIZERBASE

// Generated March 2023, David Stewart
//
//  Virtual base class used to cluster TrkrHits into TrkrClusters, but using only the reconstructed hits
//  from phg4 embedded (``truth'') tracks. In each of the following modules, a child-class will be derived
//  with the "cluster_hits()" virtual function implemented
//  - PHG4MvtxHitReco
//  - PHG4InttHitReco
//  - PHG4TpcElectronDrift
//  - (maybe?) tpot?
//
//  It's job is to:
//    (1) build TrkrTruthTracks in the TrkrTruthTrackContainer
//    (2) build TrkrClusters in the truth clusters TrkrClusterContainer
//  It does this by collecting the TrkrHit's associated with each PHG4 truth track, and when they
//  are all collected, it calls the down-stream macros (the same ones which do the Svtx clustering)
//  on this subset of the TrkrHits, and assigning these clusters to the TrkrClusterContainer and
//  the associated TrkrTruthTrackContainer.

#include <trackbase/TrkrDefs.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <map>
#include <iostream>

class PHCompositeNode;
class TrkrHitSetContainer;
class TrkrClusterContainer;
class TrkrTruthTrackContainer;
class TrkrTruthTrack;
class PHG4TruthInfoContainer;
class PHG4Hit;

class TruthClusterizerBase {
  protected:
  TrkrHitSetContainer*     m_hits          ;
  int                      m_verbosity     { 0 };
  PHCompositeNode*         m_topNode       { nullptr };
  TrkrTruthTrackContainer* m_truthtracks   { nullptr };
  TrkrClusterContainer*    m_clusters      { nullptr }; // cluster container passed to individual clusterers
  PHG4TruthInfoContainer*  m_truthinfo     { nullptr };
  int                      m_trkid         { -1      };
  bool                     m_is_emb        { false   };
  bool                     m_was_emb       { false   };
  bool                     m_is_new_track  { false   };
  TrkrTruthTrack*          m_current_track { nullptr };
  

  std::map<TrkrDefs::hitsetkey,unsigned int> m_hitsetkey_cnt {}; // counter for making ckeys form hitsetkeys

  // implemented individually for mvtx, intt and tpc cluster hits
  /* static int dummy_cluster_hits() { */
    /* return Fun4AllReturnCodes::EVENT_OK; */
  /* }; */

  public:
  TruthClusterizerBase ( );
  void    init_clusterizer_base ( PHCompositeNode*& _topNode, int verbosity );
  virtual ~TruthClusterizerBase();

  // main use functions
  void check_g4hit_status (PHG4Hit*);
  void transfer_clusters(TrkrClusterContainer*);
  void update_track();
  void transfer_clusters();

  void addhitset   (TrkrDefs::hitsetkey, TrkrDefs::hitkey, float neffelectrons);

  // convenience
  int  Verbosity() { return m_verbosity; };
  void set_verbosity(int _) { m_verbosity = _; };
  void print_clusters(int nclusprint=20);

};

#endif
