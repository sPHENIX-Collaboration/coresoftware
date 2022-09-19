// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4EVAL_TRACKEVALUATIONCONTAINERV1_H
#define G4EVAL_TRACKEVALUATIONCONTAINERV1_H

/*!
 * \file TrackEvaluationContainerv1.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TrackEvaluationContainer.h"

#include <vector>

//! track evaluation container
/*!
  Contains:
  - some global quantities, such as number of clusters in tracking detectors
  - relevant information on clusters,
  - relevant information on tracks, including redundant information about participating clusters
  Each of those can be turned on/off using flags in the filling module at runtime

  IMPORTANT NOTE: do not modify and commit neither this class or subclasses. 
  This will render past DSTs containing the container unreadible
  If you plan permanent modification to the class you need to create a TrackEvaluationContainerv2, 
  copy the relevant code from this one there and modify the TrackEvaluation module in order to use the new container
*/
class TrackEvaluationContainerv1: public TrackEvaluationContainer
{
  
  public:
  
  //! constructor
  explicit TrackEvaluationContainerv1()
  {
    // only one event structure per event (!)
    m_events.reserve(1);
  }
   
  //! reset
  void Reset() override;

  //! event information
  /*! do not modify and commit: this will break reading past DSTs */
  class EventStruct
  {

    public:
    using List = std::vector<EventStruct>;
    static constexpr size_t max_layer = 57;

    // constructor
    EventStruct()
    {
      for( size_t i = 0; i < max_layer; ++i )
      { nclusters[i] = 0; }
    }
    
    //! number of clusters per layer / event
    int nclusters[max_layer];

    //! number of clusters in the TPC
    int nclusters_mvtx = 0;

    //! number of clusters in the intt
    int nclusters_intt = 0;

    //! number of clusters in the TPC
    int nclusters_tpc = 0;

    //! number of clusters in the Micromegas
    int nclusters_micromegas = 0;
  };
  
  //! cluster information
  /*! do not modify and commit: this will break reading past DSTs */
  class ClusterStruct
  {
    public:

    using List = std::vector<ClusterStruct>;

    //! cluster layer
    unsigned int layer = 0;

    //! number of hits belonging to the cluster
    unsigned int size = 0;

    //! number of g4hits associated to this cluster
    unsigned int truth_size = 0;

    //! number of hits along phi and along z
    int phi_size = 0;
    int z_size = 0;

    //! cluster charge
    unsigned int adc = 0;

    //! cluster quality info
    int ovlp = 0;
    int edge = 0;

    //!@name cluster position
    //@{
    float x = 0;
    float y = 0;
    float z = 0;
    float r = 0;
    float phi = 0;
    float phi_error = 0;
    float z_error = 0;
    //@}

    //!@name track position at cluster
    //@{
    float trk_x = 0;
    float trk_y = 0;
    float trk_z = 0;
    float trk_r = 0;
    float trk_phi = 0;

    //! track errors
    float trk_phi_error = 0;
    float trk_z_error = 0;

    //! track inclination at cluster in r,phi plane
    float trk_alpha = 0;

    //! track inclination at cluster in r,z plane
    float trk_beta = 0;

    //! track radius
    float trk_radius = 0;

    //@}

    //!@name truth position
    //@{
    float truth_x = 0;
    float truth_y = 0;
    float truth_z = 0;
    float truth_r = 0;
    float truth_phi = 0;

    //! track inclination at cluster in r,phi plane
    float truth_alpha = 0;

    //! track inclination at cluster in r,z plane
    float truth_beta = 0;
    //@}

    //!@name charges
    //@{

    //* maximum charge on strip
    float energy_max = 0;

    //* sum of strip charges
    float energy_sum = 0;

    //@}

    //!@name track local momentum information
    //!TODO: in principle trk_alpha and trk_beta can be calculated from those. There should be no need to store them
    //@{
    float trk_px = 0;
    float trk_py = 0;
    float trk_pz = 0;
    //@}

    //!@name truth local momentum information
    //!TODO: in principle truth_alpha and truth_beta can be calculated from those. There should be no need to store them
    //@{
    float truth_px = 0;
    float truth_py = 0;
    float truth_pz = 0;
    //@}

  };

  //! track information
  /*! do not modify and commit: this will break reading past DSTs */
  class TrackStruct
  {
    public:

    // constructor
    explicit TrackStruct()
    {   
      // allocate enough size for the clusters
      static constexpr int max_layers = 60;
      clusters.reserve( max_layers );
    }

    using List = std::vector<TrackStruct>;

    int charge = 0;
    int nclusters = 0;
    int64_t mask = 0;

    int nclusters_mvtx = 0;
    int nclusters_intt = 0;
    int nclusters_tpc = 0;
    int nclusters_micromegas = 0;

    float chisquare = 0;
    int ndf = 0;

    //!@name position
    //@{
    float x = 0;
    float y = 0;
    float z = 0;
    float r = 0;
    float phi = 0;
    //@}

    //!@name momentum
    //@{
    float px = 0;
    float py = 0;
    float pz = 0;
    float pt = 0;
    float p = 0;
    float eta = 0;
    //@}

    //!@name truth momentum
    //@{
    int pid = 0;
    int embed = 0;
    bool is_primary = false;
    int gtrackID = 0;
    float truth_t = 0;

    // number of g4hits from this MC track that match
    int contributors = 0;

    float truth_px = 0;
    float truth_py = 0;
    float truth_pz = 0;
    float truth_pt = 0;
    float truth_p = 0;
    float truth_eta = 0;
    //@}

    // associate clusters
    ClusterStruct::List clusters;
  };
  
  //!@name accessors
  //@{

  const EventStruct::List& events() const 
  { return m_events; }
  
  const ClusterStruct::List& clusters() const
  { return m_clusters; }
  
  const TrackStruct::List& tracks() const
  { return m_tracks; }
  
  //@}
  
  //!@name modifiers
  //@{

  void addEvent( const EventStruct& event )
  { m_events.push_back( event ); }
  
  void addCluster( const ClusterStruct& cluster )
  { m_clusters.push_back( cluster ); }
  
  void addTrack( const TrackStruct& track )
  { m_tracks.push_back( track ); }

  void clearEvents()
  { m_events.clear(); }
  
  void clearClusters()
  { m_clusters.clear(); }
  
  void clearTracks()
  { m_tracks.clear(); }
  
  //@}
  
  private:

  //! event struct
  /* there is only one element per event in this array */
  EventStruct::List m_events;
 
  //! clusters array
  ClusterStruct::List m_clusters;
 
  //! tracks array
  TrackStruct::List m_tracks;
    
  ClassDefOverride(TrackEvaluationContainerv1,1)
    
};

#endif
