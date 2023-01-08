// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4EVAL_DSTCONTAINERTCL_H
#define G4EVAL_DSTCONTAINERTCL_H

/*!
 * \file DSTContainerTcl.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 * \modified Christof Roland <cer@mit.edu>
 */

#include "DSTContainer.h"
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrCluster.h>

#include <TObject.h>
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
  If you plan permanent modification to the class you need to create a DSTContainerv2,
  copy the relevant code from this one there and modify the DST module in order to use the new container
*/
class DSTContainerTcl: public DSTContainer
{

  public:

  //! constructor
  explicit DSTContainerTcl()
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
  class ClusterStruct : public TObject
  {
    public:

    // constructor
    explicit ClusterStruct()
    {};

    explicit ClusterStruct(TrkrDefs::cluskey key, TrkrCluster* cluster)
    {
      clusterKey = key;
      layer = TrkrDefs::getLayer(key);
      phi_seg = TrkrDefs::getPhiElement(key);
      z_seg = TrkrDefs::getZElement(key);
      loc_x = cluster->getLocalX();
      loc_y = cluster->getLocalY();
      loc_x_err = cluster->getActsLocalError(0,0);
      loc_y_err = cluster->getActsLocalError(1,1);
      adc = cluster->getAdc();

      // for v3
      int nLocal = 2;
      for (auto iLocal = 0; iLocal < nLocal; ++iLocal) {
        for (auto jLocal = 0; jLocal < nLocal; ++jLocal) {
          actsLocalError[iLocal][jLocal] =
            cluster->getActsLocalError(iLocal, jLocal);
        }
      }
      subSurfKey = cluster->getSubSurfKey();
    }

    using List = std::vector<ClusterStruct>;

    //! cluster layer
    uint8_t layer = 0;
    uint8_t phi_seg = 0;
    uint8_t z_seg = 0;
    //! number of hits belonging to the cluster
    char size = 0;

    // //! number of hits along phi and along z
    // uint8_t phi_size = 0;
    // uint8_t z_size = 0;

    //!@name cluster position
    //@{
    float loc_x = 0;
    float loc_y = 0;
    float loc_x_err = 0;
    float loc_y_err = 0;
    //@}
    // short int amp;
    //! size covariance matrix
    float cor_size[6];              //< size covariance matrix (packed storage) (+/- cm^2)
    //! error matrix
    float cor_error[6];               //< covariance matrix: rad, arc and z

    // for v3
    // float m_local[2];             //< 2D local position [cm]
    float actsLocalError[2][2];
    TrkrDefs::subsurfkey subSurfKey;
    unsigned int adc = 0;
    float rPhiError = 0;
    float zError = 0;

    // for v4
    char overlap = 0;
    char edge = 0;


    // clusterkey
    TrkrDefs::cluskey clusterKey;

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

    char charge = 0;
    char nclusters = 0;
    int64_t mask = 0;

    char nclusters_mvtx = 0;
    char nclusters_intt = 0;
    char nclusters_tpc = 0;
    char nclusters_micromegas = 0;

    float chisquare = 0;
    char ndf = 0;

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

  // get unique index in cov. matrix array from i and j
  static inline unsigned int covarIndex(unsigned int i, unsigned int j)
  {
    if (i > j) std::swap(i, j);
    return i + 1 + (j + 1) * (j) / 2 - 1;
  }

  private:

  //! event struct
  /* there is only one element per event in this array */
  EventStruct::List m_events;

  //! clusters array
  ClusterStruct::List m_clusters;

  //! tracks array
  TrackStruct::List m_tracks;

  ClassDefOverride(DSTContainerTcl,1)

};

#endif
