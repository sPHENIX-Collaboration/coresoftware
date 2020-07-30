#ifndef G4EVAL_TRACKINGEVALUATOR_HP_H
#define G4EVAL_TRACKINGEVALUATOR_HP_H

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>
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
class TrkrCluster;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;

class TrackingEvaluator_hp : public SubsysReco
{
  public:

  /// constructor
  TrackingEvaluator_hp( const std::string& = "TrackingEvaluator_hp" );

  /// global initialization
  virtual int Init(PHCompositeNode*);

  /// run initialization
  virtual int InitRun(PHCompositeNode*);

  /// event processing
  virtual int process_event(PHCompositeNode*);

  /// end of processing
  virtual int End(PHCompositeNode*);

  // event information
  class EventStruct
  {

    public:
    using List = std::vector<EventStruct>;
    static constexpr size_t max_layer = 57;

    // constructor
    EventStruct()
    {
      for( int i = 0; i < max_layer; ++i )
      {
        _nhits[i] = 0;
        _nclusters[i] = 0;
      }
    }
    int _nhits[max_layer];
    int _nclusters[max_layer];
    int _nclusters_mvtx = 0;
    int _nclusters_intt = 0;
    int _nclusters_tpc = 0;
    int _nclusters_micromegas = 0;
  };

  // cluster information to be stored in tree
  class ClusterStruct
  {
    public:

    using List = std::vector<ClusterStruct>;

    /// cluster layer
    unsigned int _layer = 0;

    /// number of hits belonging to the cluster
    unsigned int _size = 0;

    /// number of g4hits associated to this cluster
    unsigned int _truth_size = 0;

    /// number of hits along phi and along z
    int _phi_size = 0;
    int _z_size = 0;

    ///@name cluster position
    //@{
    float _x = 0;
    float _y = 0;
    float _z = 0;
    float _r = 0;
    float _phi = 0;
    float _phi_error = 0;
    float _z_error = 0;
    //@}

    ///@name track position at cluster
    //@{
    float _trk_x = 0;
    float _trk_y = 0;
    float _trk_z = 0;
    float _trk_r = 0;
    float _trk_phi = 0;

    /// track errors
    float _trk_phi_error = 0;
    float _trk_z_error = 0;

    /// track inclination at cluster in r,phi plane
    float _trk_alpha = 0;

    /// track inclination at cluster in r,z plane
    float _trk_beta = 0;

    //@}

    ///@name truth position
    //@{
    float _truth_x = 0;
    float _truth_y = 0;
    float _truth_z = 0;
    float _truth_r = 0;
    float _truth_phi = 0;

    /// track inclination at cluster in r,phi plane
    float _truth_alpha = 0;

    /// track inclination at cluster in r,z plane
    float _truth_beta = 0;
    //@}

  };

  // track information to be stored in tree
  class TrackStruct
  {
    public:

    // constructor
    explicit TrackStruct()
    {
      // allocate enough size for the clusters
      static constexpr int max_layers = 60;
      _clusters.reserve( max_layers );
    }

    using List = std::vector<TrackStruct>;

    int _charge = 0;
    int _nclusters = 0;
    int64_t _mask = 0;

    int _nclusters_mvtx = 0;
    int _nclusters_intt = 0;
    int _nclusters_tpc = 0;
    int _nclusters_micromegas = 0;

    float _chisquare = 0;
    int _ndf = 0;

    ///@name position
    //@{
    float _x = 0;
    float _y = 0;
    float _z = 0;
    float _r = 0;
    float _phi = 0;
    //@}

    ///@name momentum
    //@{
    float _px = 0;
    float _py = 0;
    float _pz = 0;
    float _pt = 0;
    float _p = 0;
    float _eta = 0;
    //@}

    ///@name truth momentum
    //@{
    int _pid = 0;
    int _embed = 0;
    bool _is_primary = false;

    // number of g4hits from this MC track that match
    int _contributors = 0;

    float _truth_px = 0;
    float _truth_py = 0;
    float _truth_pz = 0;
    float _truth_pt = 0;
    float _truth_p = 0;
    float _truth_eta = 0;
    //@}

    // associate clusters
    ClusterStruct::List _clusters;
  };

  // pair information to be stored in tree
  class TrackPairStruct
  {
    public:

    using List = std::vector<TrackPairStruct>;

    int _charge = 0;

    ///@name momentum
    //@{
    float _px = 0;
    float _py = 0;
    float _pz = 0;
    float _pt = 0;
    float _p = 0;
    float _e = 0;
    float _m = 0;
    float _eta = 0;

    std::array<float,2> _trk_pt = {{0,0}};

    //@}
  };

  /// track container
  class Container: public PHObject
  {

    public:

    /// constructor
    explicit Container() = default;

    /// copy constructor
    explicit Container(const Container &) = delete;

    /// assignment operator
    Container& operator = ( const Container& ) = delete;

    /// reset
    virtual void Reset();

    ///@name accessors
    //@{

    const EventStruct::List& events() const
    { return _events; }

    const ClusterStruct::List& clusters() const
    { return _clusters; }

    const TrackStruct::List& tracks() const
    { return _tracks; }

    const TrackPairStruct::List& trackPairs() const
    { return _track_pairs; }

    //@}

    ///@name modifiers
    //@{

    void addEvent( const EventStruct& event )
    { _events.push_back( event ); }

    void addCluster( const ClusterStruct& cluster )
    { _clusters.push_back( cluster ); }

    void addTrack( const TrackStruct& track )
    { _tracks.push_back( track ); }

    void addTrackPair( const TrackPairStruct& pair )
    { _track_pairs.push_back( pair ); }

    void clearEvents()
    { _events.clear(); }

    void clearClusters()
    { _clusters.clear(); }

    void clearTracks()
    { _tracks.clear(); }

    void clearTrackPairs()
    { _track_pairs.clear(); }

    //@}

    private:

    /// event struct
    EventStruct::List _events;

    /// clusters array
    ClusterStruct::List _clusters;

    /// tracks array
    TrackStruct::List _tracks;

    /// track pairs array
    TrackPairStruct::List _track_pairs;

    ClassDef(Container,1)

  };

  enum Flags
  {

    EvalEvent = 1<<0,
    EvalClusters = 1<<1,
    EvalTracks = 1<<2,
    EvalTrackPairs = 1<<3,
    PrintClusters = 1<<4,
    PrintTracks = 1<<5
  };

  /// set flags. Should be a bitwise or of Flags enum
  void set_flags( int flags )
  { m_flags = flags; }

  private:

  /// load nodes
  int load_nodes( PHCompositeNode* );

  /// evaluate evens
  void evaluate_event();

  /// evaluate clusters
  void evaluate_clusters();

  /// evaluate tracks
  void evaluate_tracks();

  /// evaluate track pairs
  void evaluate_track_pairs();

  /// print clusters
  void print_clusters() const;

  /// print tracks
  void print_tracks() const;

  /// print cluster and association
  void print_cluster( TrkrDefs::cluskey, TrkrCluster* ) const;

  /// print track content
  void print_track( SvtxTrack* ) const;

  // get geant hits associated to a cluster
  using G4HitSet = std::set<PHG4Hit*>;
  G4HitSet find_g4hits( TrkrDefs::cluskey ) const;

  // get G4Particle id of max contributor to a given track
  std::pair<int,int> get_max_contributor( SvtxTrack* ) const;

  // get embedded id for given g4track
  int get_embed(PHG4Particle*) const;

  // cluster array
  Container* m_container = nullptr;

  // flags
  int m_flags = EvalEvent | EvalClusters | EvalTracks | EvalTrackPairs;

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

};

#endif  // G4EVAL_TRACKINGEVALUATOR_HP_H
