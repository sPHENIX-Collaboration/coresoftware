// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4EVAL_DSTCONTAINERV3_H
#define G4EVAL_DSTCONTAINERV3_H

/*!
 * \file DSTContainerv3.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 * \modified Christof Roland <cer@mit.edu>
 */

#include "DSTContainer.h"
#include <TObject.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterv4.h>

#include <vector>
#include <TClonesArray.h>

// #include <boost/smart_ptr.hpp>

class SvtxTrack;
class TrkrCluster;
class TrkrClusterv4;

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
class DSTContainerv3: public DSTContainer
{

  public:

  //! constructor
  explicit DSTContainerv3()
  {
    // only one event structure per event (!)
    m_events.reserve(1);
    std::cout << "init array" << "\n";
    trkrClsDST = new TClonesArray("TrkrClusterv4");
    arrClsDST = new TClonesArray("DSTContainerv3::ClusterStruct");
    arrTrkDST = new TClonesArray("DSTContainerv3::TrackStruct");
    // arrKeyDST = new TClonesArray("TrkrDefs::hitsetkey");
    arrKeyDST = new TClonesArray("DSTContainerv3::ClusterKeyStruct");
    std::cout << " end init array" << std::endl;
  }

  //! reset
  void Reset() override;

  //! event information
  /*! do not modify and commit: this will break reading past DSTs */
  class EventStruct : public PHObject
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
  class ClusterStruct : public PHObject
  {
    public:

    // constructor
    ClusterStruct() {};

    // ClusterStruct(TrkrDefs::hitsetkey hitsetkeyIn, TrkrClusterv4 clusIn)
    ClusterStruct(TrkrDefs::hitsetkey hitsetkeyIn, TrkrClusterv4 clusIn, TrkrDefs::cluskey ckeyIn)
    // ClusterStruct(TrkrDefs::hitsetkey hitsetkeyIn)
    {
      clusterKey = ckeyIn;
      hitsetkey = hitsetkeyIn;
      cluster = clusIn;
    }

    ClusterStruct(TrkrDefs::hitsetkey hitsetkeyIn, TrkrClusterv4 clusIn)
    {
      hitsetkey = hitsetkeyIn;
      cluster = clusIn;
    }

    using List = std::vector<ClusterStruct>;


    // clusterkey
    TrkrDefs::cluskey clusterKey;
    TrkrDefs::hitsetkey hitsetkey;
    TrkrClusterv4 cluster;

  };

  class ClusterKeyStruct : public PHObject
  {
  public:

    // constructor
    ClusterKeyStruct() {};

    ClusterKeyStruct(TrkrDefs::hitsetkey hitsetkeyIn)
    {
      hitsetkey = hitsetkeyIn;
    }

    // clusterkey
    TrkrDefs::hitsetkey hitsetkey;
  };

  //! track information
  /*! do not modify and commit: this will break reading past DSTs */
  class TrackStruct : public PHObject
  {
    public:

    // constructor
    explicit TrackStruct()
    {
      // allocate enough size for the clusters
      static constexpr int max_layers = 60;
      clusters.reserve( max_layers );
    }

    // TrackStruct(SvtxTrack track)
    // {
    //   // allocate enough size for the clusters
    //   static constexpr int max_layers = 60;
    //   clusters.reserve(max_layers);
    //   set_id(track.get_id());
    // }

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

  // typedef boost::shared_ptr<TClonesArray> arr_ptr;
  TClonesArray* arrEvtDST = nullptr;
  TClonesArray* arrTrkDST = nullptr;
  TClonesArray* arrClsDST = nullptr;
  TClonesArray* trkrClsDST = nullptr;
  TClonesArray* arrKeyDST = nullptr;
  // TClonesArray* arrClsDST = new TClonesArray("DSTContainerv3::ClusterStruct");

  private:

  //! event struct
  /* there is only one element per event in this array */
  EventStruct::List m_events;

  //! clusters array
  ClusterStruct::List m_clusters;

  //! tracks array
  TrackStruct::List m_tracks;

  // expt
  // std::vector<TrkrClusterv4*> clus;
  // TrkrClusterv4* clus;

  ClassDefOverride(DSTContainerv3,1)

};

#endif
