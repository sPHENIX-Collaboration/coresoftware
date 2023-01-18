#ifndef G4TPC_TPCCLUSTERBUILDER_H
#define G4TPC_TPCCLUSTERBUILDER_H

// Revised 04-Dec-2022, David Stewart
// basic use: 
//  (a) optionally used set_current_track() (pointless to not set the track though,
//      otherwise the new TrkrClsuters won't be assigned to the proper location)
//  (b) fill with TrkrHits in addhitset() (in the MapToPadPlane module)
//  (c) after all TrkrHit's have been added for a given truth track, 
//      call cluster_and_reset(false) to generate the TrkrClusters from the
//      TrkrHit's, fill in the clusterkeys in the truth_track, and reset
//      the TrkrHit container.
//  (d) at the end of the event, run cluster_and_reset(true) which will check
//      if there are clusters to make, and also clear out the counter for the
//      hitsetkeys
//
//  Note:
//  - the algorithm to get the local surface will probably be updated in 
//    TpcClusterizer, when that happens update here as well

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/ActsGeometry.h>
#include <map>
#include <climits>

class TrkrCluster;
class PHG4TpcCylinderGeom;
class TrkrClusterContainer;
class TrkrHitSetContainer;
class TrkrTruthTrack;
class PHG4TpcCylinderGeomContainer;
class TrkrTruthTrackContainer;

// This is the basic data for each set of TrkrHits from each TrkrHitsSet 
// to be used in tpc/TpcClusterizer.cc
class TpcClusterBuilder {
  double square(double);
  double square(float);

  TrkrClusterContainer*         m_clusterlist; // fill for output
  ActsGeometry*                 m_tGeometry;   // used to generate clusters
  PHG4TpcCylinderGeomContainer* geom_container;

  // internal containers to fill and consume hits and fill with tracks
  TrkrHitSetContainer* m_hits        { new TrkrHitSetContainerv1() };
  TrkrTruthTrack*      current_track { nullptr };
  std::map<TrkrDefs::hitsetkey,unsigned int> hitsetkey_cnt {};

  int n_tracks {0};

  int verbosity {0};

  public:
  TpcClusterBuilder( 
      TrkrClusterContainer*         _truth_cluster_container
    , ActsGeometry*                 _ActsGeometry
    , PHG4TpcCylinderGeomContainer* _geom_container
  );

  bool is_embedded_track {false};
  void cluster_and_reset (bool clear_hitsetkey_cnt);
  void addhitset (TrkrDefs::hitsetkey, TrkrDefs::hitkey, float neffelectrons);
  void set_current_track (TrkrTruthTrack* _trkrtruthtrack);
  void print(TrkrTruthTrackContainer*, int nclusprint=-1);
  void print_file(TrkrTruthTrackContainer*, std::string);
  void set_verbosity(int verbosity_level);

  ~TpcClusterBuilder(){
    delete m_hits;
  };

  private: // from TpcClusterizer.h parameters; also used as general
  /* bool   do_wedge_emulation = false; */
  /* double SectorFiducialCut = 0.5; */
  /* unsigned short MaxClusterHalfSizePhi = 3; */
  /* unsigned short MaxClusterHalfSizeT   = 5; */
  /* int    cluster_version = 4; */
  double AdcClockPeriod  = 53.0; // ns

  // TPC shaping offset correction parameter
  // From Tony Frawley July 5, 2022
  /* double m_sampa_tbias = 39.6;  // ns */  

  void reset(bool clear_hitsetkey_cnt);
};

#endif  //TRACKBASE_PADPLANEREADOUTSTRUCT_H
