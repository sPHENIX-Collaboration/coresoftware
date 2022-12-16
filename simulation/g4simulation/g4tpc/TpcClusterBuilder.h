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

  TrkrClusterContainer*         m_clusterlist; // fill for output
  ActsGeometry*                 m_tGeometry;               // used to generate clusters
  PHG4TpcCylinderGeomContainer* geom_container;

  // internal containers to fill and consume hits and fill with tracks
  TrkrHitSetContainer* m_hits { new     TrkrHitSetContainerv1() };
  TrkrTruthTrack*      current_track   { nullptr };
  std::map<TrkrDefs::hitsetkey,unsigned int> hitsetkey_cnt {};

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
  void cluster_and_reset (bool clear_hitsetkey_cnt);
  void addhitset (TrkrDefs::hitsetkey, TrkrDefs::hitkey, float neffelectrons);
  void set_current_track (TrkrTruthTrack* _trkrtruthtrack);


  private: // from TpcClusterizer.h parameters; also used as general
  /* bool   do_wedge_emulation = false; */
  /* double pedestal          = 74.4; */
  double SectorFiducialCut = 0.5;
  /* unsigned short MaxClusterHalfSizePhi = 3; */
  /* unsigned short MaxClusterHalfSizeT   = 5; */
  /* int    cluster_version = 4; */
  double AdcClockPeriod  = 53.0; // ns

  // TPC shaping offset correction parameter
  // From Tony Frawley July 5, 2022
  double m_sampa_tbias = 39.6;  // ns  

  ~TpcClusterBuilder(){
    delete m_hits;
  };
};

#endif  //TRACKBASE_PADPLANEREADOUTSTRUCT_H
