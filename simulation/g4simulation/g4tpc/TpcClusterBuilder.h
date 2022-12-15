#ifndef G4TPC_TPCCLUSTERBUILDER_H
#define G4TPC_TPCCLUSTERBUILDER_H

// Revised 04-Dec-2022, David Stewart
//
//
// // put in a switch for skipping the noise subtraction
//
// Essentially a local clone of TpcClusterizer.{cc,h} into the wrapper object TpcClusterbuilder,
// which object use is as follows:
//
//  Initialize with: 
//    - output TrkrClusterContainer, where the generated clusters will be put
//    - ActsGeometry and PHG4TpcCylinderGeomContainer
//  Optionally update with current TrkrTruthTrack which will collect TrkrDefs::cluskey from
//     the TrkrCluster's put in the container
//
// basic use: 
//  (a) optionally used set_current_track()
//  (b) fill with TrkrHits in addhitset() (in the MapToPadPlane module)
//  (c) after all TrkrHit's have been added, call cluster_and_reset()
//      to generate the TrkrClusters from the TrkrHit's and reset the container
//
// future use/room for improvement:
//  - the posix threading from TpcClsuterizer.{cc,h} have been removed. For first used 
//    case this is probably ok, as it is run on a single tracks set of hits at a time
//    which is a relatively small set
//  - there is majority code overlap in literally copied code TpcClusterizer.{cc,h} which
//    immediately suggests a few things:
//      (a) extracting the module out to a common location (perhaps in TpcClusterizer)
//          so that there will be only one copy of code
//      (b) as an alternative to (a), carefully simplify this local code
//      (c) add a switch to write out the multimap as well
//

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/ActsGeometry.h>
#include <map>
#include <climits>

class TrkrCluster;
class PHG4TpcCylinderGeom;
/* class TrkrTruthTrackContainer; */
class TrkrClusterContainer;
class TrkrHitSetContainer;
class TrkrTruthTrack;
class PHG4TpcCylinderGeomContainer;

  // This is really just a structure used in simulation/g4simulation/g4tpc/PHG4TpcElectronDrift.cc
  // to collect necessary statistics from simulation/g4simulation/g4tpc/PHG4TpcPadPlaneReadout.cc
  // 
  // It could have been just std::pair< std::array<int,6>, std::pair<double,double>>, but 
  // having the structure to name those subtypes.

// This is the basic data for each set of TrkrHits from each TrkrHitsSet 
// to be used in tpc/TpcClusterizer.cc
class TpcClusterBuilder {
  double square(double);
  double square(float);

  struct thread_data 
  {
    PHG4TpcCylinderGeom *layergeom       = nullptr;
    TrkrHitSet          *hitset          = nullptr;
    /* RawHitSetv1         *rawhitset       = nullptr; */
    ActsGeometry   *tGeometry          = nullptr;
    unsigned int    layer              = 0;
    int             side               = 0;
    unsigned int    sector             = 0;
    float           pedestal           = 0;
    bool            do_wedge_emulation = true;
    unsigned short  phibins            = 0;
    unsigned short  phioffset          = 0;
    unsigned short  tbins              = 0;
    unsigned short  toffset            = 0;
    unsigned short  maxHalfSizeT       = 0;
    unsigned short  maxHalfSizePhi     = 0;
    double          m_tdriftmax        = 0;
    double          sampa_tbias        = 0;
    int             cluster_version    = 4;
    /* std::vector<assoc>   association_vector; */
    std::vector<TrkrCluster*> cluster_vector;
    int                  verbosity       = 0;

    bool            m_skip_noise       = false;
  };

  // basic structure used to pass data in tpc/TpcClusterizer.cc
  struct ihit
  {
    unsigned short iphi = 0;
    unsigned short it = 0;
    unsigned short adc = 0;
    unsigned short edge = 0;
  };

  // pointers to objects passed from the outside
  TrkrClusterContainer*         m_clusterlist; // fill for output
  /* TrkrTruthTrackContainer*      truth_track_container;   // fill for output */
  ActsGeometry*                 m_tGeometry;               // used to generate clusters
  PHG4TpcCylinderGeomContainer* geom_container;

  // internal containers to fill and consume hits and fill with tracks
  TrkrHitSetContainer* m_hits { new     TrkrHitSetContainerv1() };
  TrkrTruthTrack*      current_track   { nullptr };

  public:
  TpcClusterBuilder( 
      TrkrClusterContainer*         _truth_cluster_container
    /* , TrkrTruthTrackContainer*      _truth_track_container */
    , ActsGeometry*                 _ActsGeometry
    , PHG4TpcCylinderGeomContainer* _geom_container
  );

  // basic use: 
  //  (a) set_current_track()
  //  (b) fill with TrkrHits in addhitset() (in the MapToPadPlane module)
  //  (c) after all TrkrHit's have been added, call cluster_and_reset()
  //      to generate the TrkrClusters from the TrkrHit's and reset the container
  //  Use is_embedded_track to switch when TrkrHit's are or are-not collected

  bool is_embedded_track {false};
  void cluster_and_reset ();
  void addhitset (TrkrDefs::hitsetkey, TrkrDefs::hitkey, float neffelectrons);
  void set_current_track (TrkrTruthTrack* _trkrtruthtrack);


  private: // from TpcClusterizer.h parameters; also used as general
  // internal functions, called inside of cluster_and_reset
  void ProcessSectorData(thread_data* my_data);
  void calc_cluster_parameter( const std::vector<ihit> &ihit_list, 
      thread_data& my_data, int ntouch, int nedge );
  void remove_hits(
      std::vector<ihit> &ihit_list, 
      std::multimap<unsigned short, ihit> &all_hit_map,
      std::vector<std::vector<unsigned short>> &adcval);
  void remove_hit(double adc, int phibin, int tbin, int edge,
      std::multimap<unsigned short, ihit> &all_hit_map,
      std::vector<std::vector<unsigned short>> &adcval);
  void get_cluster(int phibin, int tbin, 
      const thread_data& my_data, const std::vector<std::vector<unsigned short>> &adcval, 
      std::vector<ihit> &ihit_list, int &touch, int &edge);
  void find_t_range(int phibin, int tbin, const thread_data& my_data, const
      std::vector<std::vector<unsigned short>> &adcval, int& tdown, int& tup,
      int &touch, int &edge);
  void find_phi_range(int phibin, int tbin, const thread_data& my_data, const
      std::vector<std::vector<unsigned short>> &adcval, int& phidown, int&
      phiup, int &touch, int &edge);

  /* TrkrHitSetContainer  *m_hits            = nullptr; */
  /* RawHitSetContainer   *m_rawhits         = nullptr; */
  /* TrkrClusterContainer *m_clusterlist     = nullptr; */
  /* TrkrClusterHitAssoc  *m_clusterhitassoc = nullptr; */
  /* ActsGeometry         *m_tGeometry       = nullptr; */
  /* bool   do_hit_assoc       = true; */
  bool   do_wedge_emulation = false;
  /* bool   do_sequential      = false; */
  /* bool   do_read_raw        = false; */
  /* bool   do_singles         = false; */
  double pedestal          = 74.4;
  double SectorFiducialCut = 0.5;
  unsigned short MaxClusterHalfSizePhi = 3;
  unsigned short MaxClusterHalfSizeT   = 5;
  int    cluster_version = 4;
  double m_tdriftmax     = 0;
  double AdcClockPeriod  = 53.0; // ns

  // TPC shaping offset correction parameter
  // From Tony Frawley July 5, 2022
  double m_sampa_tbias = 39.6;  // ns  

  ~TpcClusterBuilder(){
    delete m_hits;
  };
};

#endif  //TRACKBASE_PADPLANEREADOUTSTRUCT_H
