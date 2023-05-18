#ifndef G4EVAL_TRACKINGEVALUATOR_HP_H
#define G4EVAL_TRACKINGEVALUATOR_HP_H

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>
#include <tpc/TpcDistortionCorrection.h>
#include <tpc/TpcClusterZCrossingCorrection.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/ActsTransformations.h>

#include <map>
#include <set>
#include <string>
#include <vector>

class ActsGeometry;
class CMFlashClusterContainer;
class PHG4CylinderGeomContainer;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class PHG4TpcCylinderGeomContainer;
class PHG4TruthInfoContainer;
class SvtxTrack;
class SvtxTrackMap;
class SvtxTrackState;
class TrkrCluster;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;

class TrackingEvaluator_hp : public SubsysReco
{
  public:

  //! constructor
  TrackingEvaluator_hp( const std::string& = "TrackingEvaluator_hp" );

  //! global initialization
  virtual int Init(PHCompositeNode*);

  //! run initialization
  virtual int InitRun(PHCompositeNode*);

  //! event processing
  virtual int process_event(PHCompositeNode*);

  //! end of processing
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
      for( size_t i = 0; i < max_layer; ++i )
      {
        _nhits[i] = 0;
        _nhits_raw[i] = 0;
        _nclusters[i] = 0;

      }
    }

    //! number of hits per layer / event
    /* for TPC hits, the charge is compared to threshold */
    unsigned int _nhits[max_layer];

    //! number of hits per layer / event
    unsigned int _nhits_raw[max_layer];

    //! number of clusters per layer / event
    unsigned int _nclusters[max_layer];

    //! number of clusters in the mvtx
    unsigned int _nclusters_mvtx = 0;

    //! number of clusters in the intt
    unsigned int _nclusters_intt = 0;

    //! number of clusters in the TPC
    unsigned int _nclusters_tpc = 0;

    //! number of clusters in the Micromegas
    unsigned int _nclusters_micromegas = 0;
    
    //! number of g4hits in the mvtx
    unsigned int _ng4hits_mvtx = 0;

    //! number of g4hits in the intt
    unsigned int _ng4hits_intt = 0;

    //! number of g4hits in the TPC
    unsigned int _ng4hits_tpc = 0;

    //! number of g4hits in the Micromegas
    unsigned int _ng4hits_micromegas = 0;
   
    //! number of central membrane clusters
    unsigned int _ncmclusters = 0;
    
  };

  // cluster information to be stored in tree
  class ClusterStruct
  {
    public:

    using List = std::vector<ClusterStruct>;

    //! cluster layer
    unsigned int _layer = 0;

    //! number of hits belonging to the cluster
    unsigned int _size = 0;

    //! number of g4hits associated to this cluster
    unsigned int _truth_size = 0;

    //! number of particles associated to this cluster
    unsigned int _ncontributors = 0;


    //! number of hits along phi and along z
    unsigned int _phi_size = 0;
    unsigned int _z_size = 0;

    //!@name cluster position
    //@{
    float _x = 0;
    float _y = 0;
    float _z = 0;
    float _r = 0;
    float _phi = 0;
    float _phi_error = 0;
    float _z_error = 0;
    //@}

    //!@name track position at cluster
    //@{
    float _trk_x = 0;
    float _trk_y = 0;
    float _trk_z = 0;
    float _trk_r = 0;
    float _trk_phi = 0;

    //! track errors
    float _trk_phi_error = 0;
    float _trk_z_error = 0;

    // extrapolation distance
    float _trk_dr = 0;
    
    //! track inclination at cluster in r,phi plane
    float _trk_alpha = 0;

    //! track inclination at cluster in r,z plane
    float _trk_beta = 0;

    //@}

    //!@name truth position
    //@{
    float _truth_x = 0;
    float _truth_y = 0;
    float _truth_z = 0;
    float _truth_r = 0;
    float _truth_phi = 0;

    //! track inclination at cluster in r,phi plane
    float _truth_alpha = 0;

    //! track inclination at cluster in r,z plane
    float _truth_beta = 0;
    //@}

    //!@name charges
    //@{

    //* maximum charge on strip
    float _energy_max = 0;

    //* sum of strip charges
    float _energy_sum = 0;

    //@}

    //!@name track local momentum information
    //!TODO: in principle trk_alpha and trk_beta can be calculated from those. There should be no need to store them
    //@{
    float _trk_px = 0;
    float _trk_py = 0;
    float _trk_pz = 0;
    //@}

    //!@name truth local momentum information
    //!TODO: in principle truth_alpha and truth_beta can be calculated from those. There should be no need to store them
    //@{
    float _truth_px = 0;
    float _truth_py = 0;
    float _truth_pz = 0;
    //@}

    //! micromegas tile id
    int _tileid = 0;

  };


  // cluster information to be stored in tree
  class CMClusterStruct
  {
    public:

    using List = std::vector<CMClusterStruct>;

    //! number of participating clusters
    unsigned int _nclusters = 0;

    //!@name cluster position
    //@{
    float _x = 0;
    float _y = 0;
    float _z = 0;
    float _r = 0;
    float _phi = 0;
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
    unsigned int _nclusters = 0;

    /// mask of layers for which there is a cluster in the track
    int64_t _mask = 0LL;

    /// mask of layers for which there is a cluster in the track, with correct associated g4hit
    int64_t _correct_mask = 0LL;

    /// mask of layers for which there is a cluster in the track, with correct associated g4hit
    int64_t _correct_mask_strict = 0LL;

    /// maks of layers for which there is a g4hit in the track
    int64_t _truth_mask = 0LL;

    unsigned int _nclusters_mvtx = 0;
    unsigned int _nclusters_intt = 0;
    unsigned int _nclusters_tpc = 0;
    unsigned int _nclusters_micromegas = 0;

    float _chisquare = 0;
    int _ndf = 0;

    //!@name position
    //@{
    float _x = 0;
    float _y = 0;
    float _z = 0;
    float _r = 0;
    float _phi = 0;
    //@}

    //!@name momentum
    //@{
    float _px = 0;
    float _py = 0;
    float _pz = 0;
    float _pt = 0;
    float _p = 0;
    float _eta = 0;
    //@}

    //!@name truth momentum
    //@{
    int _pid = 0;
    int _embed = 0;
    bool _is_primary = false;

    // number of g4hits from this MC track that match
    unsigned int _contributors = 0;

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
  
  
  // track information to be stored in tree
  class TrackStruct_small
  {
    public:

    // constructor
    explicit TrackStruct_small() = default;
 
    int _charge = 0;
    unsigned int _nclusters = 0;

    /// mask of layers for which there is a cluster in the track
    int64_t _mask = 0LL;

    unsigned int _nclusters_mvtx = 0;
    unsigned int _nclusters_intt = 0;
    unsigned int _nclusters_tpc = 0;
    unsigned int _nclusters_micromegas = 0;

    float _chisquare = 0;
    int _ndf = 0;

    float _px = 0;
    float _py = 0;
    float _pz = 0;
    float _pt = 0;
    float _p = 0;
    float _eta = 0;

  };

  // pair information to be stored in tree
  class TrackPairStruct
  {
    public:

    using List = std::vector<TrackPairStruct>;

    int _charge = 0;

    //!@name momentum
    //@{
    float _px = 0;
    float _py = 0;
    float _pz = 0;
    float _pt = 0;
    float _p = 0;
    float _e = 0;
    float _m = 0;
    float _eta = 0;

    std::array<TrackStruct_small,2> _tracks;
    
    //@}
  };

  //! track container
  class Container: public PHObject
  {

    public:

    //! constructor
    explicit Container() = default;

    //! copy constructor
    explicit Container(const Container &) = delete;

    //! assignment operator
    Container& operator = ( const Container& ) = delete;

    //! reset
    virtual void Reset();

    //!@name accessors
    //@{

    const EventStruct::List& events() const
    { return _events; }

    const ClusterStruct::List& clusters() const
    { return _clusters; }

    const CMClusterStruct::List& cm_clusters() const 
    { return _cm_clusters; }
    
    const TrackStruct::List& tracks() const
    { return _tracks; }

    const TrackPairStruct::List& trackPairs() const
    { return _track_pairs; }

    //@}

    //!@name modifiers
    //@{

    void addEvent( const EventStruct& event )
    { _events.push_back( event ); }

    void addCluster( const ClusterStruct& cluster )
    { _clusters.push_back( cluster ); }

    void addCMCluster( const CMClusterStruct& cluster )
    { _cm_clusters.push_back( cluster ); }
    
    void addTrack( const TrackStruct& track )
    { _tracks.push_back( track ); }

    void addTrackPair( const TrackPairStruct& pair )
    { _track_pairs.push_back( pair ); }

    void clearEvents()
    { _events.clear(); }

    void clearClusters()
    { _clusters.clear(); }

    void clearCMClusters()
    { _cm_clusters.clear(); }

    void clearTracks()
    { _tracks.clear(); }

    void clearTrackPairs()
    { _track_pairs.clear(); }

    //@}

    private:

    //! event struct
    EventStruct::List _events;

    //! clusters array
    ClusterStruct::List _clusters;

    //! central membrane clustsrs
    CMClusterStruct::List _cm_clusters;
    
    //! tracks array
    TrackStruct::List _tracks;

    //! track pairs array
    TrackPairStruct::List _track_pairs;

    ClassDef(Container,1)

  };

  enum Flags
  {
    EvalEvent = 1<<0,
    EvalClusters = 1<<1,
    EvalCMClusters = 1<<6,
    EvalTracks = 1<<2,
    EvalTrackPairs = 1<<3,
    PrintClusters = 1<<4,
    PrintTracks = 1<<5
  };

  //! set flags. Should be a bitwise or of Flags enum
  void set_flags( int flags )
  { m_flags = flags; }

  //! tracl map name
  void set_trackmapname( const std::string& value )
  { m_trackmapname = value; }
  
  //! utility functions
  static bool has_layer( int64_t mask, int layer )
  { return mask & (1LL<<layer); }

  //! get number of clusters in given range
  static int get_nclusters( int64_t mask, int first, int last )
  {
    int out = 0;
    for( int layer = first; layer < last; ++layer )
    { out += (int) has_layer( mask, layer ); }

    return out;
  }

  //! get number of mvtx clusters from mask
  static int get_nclusters_mvtx( int64_t mask )
  { return get_nclusters( mask, 0, 3 ); }

  //! get number of intt clusters from mask
  static int get_nclusters_intt( int64_t mask )
  { return get_nclusters( mask, 3, 7 ); }

  //! get number of tpc clusters from mask
  static int get_nclusters_tpc( int64_t mask )
  { return get_nclusters( mask, 7, 55 ); }

  //! get number of micromegas clusters from mask
  static int get_nclusters_micromegas( int64_t mask )
  { return get_nclusters( mask, 55, 57 ); }

  /// cluster version
  /* Note: this could be retrived automatically using dynamic casts from TrkrCluster objects */
  void set_cluster_version(int value) { m_cluster_version = value; }

  private:

  //! load nodes
  int load_nodes( PHCompositeNode* );

  //! evaluate event
  void evaluate_event();

  //! evaluate clusters
  void evaluate_clusters();

  //! evaluate cm clusters
  void evaluate_cm_clusters();
  
  //! evaluate tracks
  void evaluate_tracks();

  //! evaluate track pairs
  void evaluate_track_pairs();

  //! print clusters
  void print_clusters() const;

  //! print tracks
  void print_tracks() const;

  //! print cluster and association
  void print_cluster( TrkrDefs::cluskey, TrkrCluster* ) const;

  //! print track content
  void print_track( SvtxTrack* ) const;

  // get geant hits associated to a cluster
  using G4HitSet = std::set<PHG4Hit*>;
  G4HitSet find_g4hits( TrkrDefs::cluskey ) const;

  //! get G4Particle id of max contributor to a given track
  std::pair<int,int> get_max_contributor( SvtxTrack* ) const;

  //! get embedded id for given g4track
  int get_embed(PHG4Particle*) const;

  //! create cluster structure from cluster
  ClusterStruct create_cluster( TrkrDefs::cluskey, TrkrCluster* ) const;

  //! add track information to a cluster
  void add_trk_information( ClusterStruct&, SvtxTrackState* ) const;

  //! add track information to a cluster for the micromegas case
  /*!
   * the difference between this and the generic method is that the track state to
   * the tiles detector plane, and not to the same radius as the cluster
   */
  void add_trk_information_micromegas( ClusterStruct&, SvtxTrackState* ) const;

  // add truth information
  void add_truth_information( ClusterStruct&, std::set<PHG4Hit*> ) const;

  // add truth information
  /*!
   * the difference between this and the generic method is that the track state to
   * the tiles detector plane, and not to the same radius as the cluster
   */
  void add_truth_information_micromegas( ClusterStruct&, std::set<PHG4Hit*> ) const;

  //! fill MC track map
  void fill_g4particle_map();

  /// get global position for a given cluster
  /**
   * uses ActsTransformation to convert cluster local position into global coordinates
   */
  Acts::Vector3 get_global_position(TrkrDefs::cluskey, TrkrCluster*, short int crossing = 0) const;

  //! evaluation node
  Container* m_container = nullptr;

  //! flags
  int m_flags = EvalEvent | EvalClusters | EvalTracks | EvalTrackPairs;

  /// Acts tracking geometry for surface lookup
  ActsGeometry *m_tGeometry = nullptr;

  // crossing z correction
  TpcClusterZCrossingCorrection m_clusterCrossingCorrection;
  
  // distortion corrections
  TpcDistortionCorrectionContainer* m_dcc_static = nullptr;
  TpcDistortionCorrectionContainer* m_dcc_average = nullptr;
  TpcDistortionCorrectionContainer* m_dcc_fluctuation = nullptr;

  /// tpc distortion correction utility class
  TpcDistortionCorrection m_distortionCorrection;
  
  //! hits
  TrkrHitSetContainer* m_hitsetcontainer = nullptr;

  //! clusters
  TrkrClusterContainer* m_cluster_map = nullptr;

  //! central membrane clusters
  CMFlashClusterContainer* m_cm_cluster_map = nullptr;

  //! cluster to hit association
  TrkrClusterHitAssoc* m_cluster_hit_map = nullptr;

  //! hit to truth association
  TrkrHitTruthAssoc* m_hit_truth_map = nullptr;

  //! track map name
  std::string m_trackmapname = "SvtxTrackMap";
  
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
  
  //! map cluster keys to g4hits
  using G4HitMap = std::map<TrkrDefs::cluskey,G4HitSet>;
  mutable G4HitMap m_g4hit_map;

  //! map trk_id to layer mask
  /** copied from SimEvaluator_hp */
  using G4ParticleMap = std::map<int,int64_t>;
  G4ParticleMap m_g4particle_map;
  
  /// cluster error parametrisation
  ClusterErrorPara m_cluster_error_parametrization;
  
  /// cluster version
  int m_cluster_version = 4;
  
};

#endif  // G4EVAL_TRACKINGEVALUATOR_HP_H
