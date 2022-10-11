#ifndef TRUTHTRKMATCHER__H
#define TRUTHTRKMATCHER__H

#include <trackbase/TrkrDefs.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

// needed:
// Nodes: 
//  Get TrkrTruthTrack's node
//  Get reconstructed tracks node
//  Compare TrkrTruthTracks::phi,eta,(pT?) vs reconstructed tracks::phi,eta,(pT?)
//    -- if there are multiple matches, then go and match the TrkrClusters (from the TrkrClusterContainer)
//  Write values to the MatchTracks on the tree:
//    nClusterHits, nClusterHitsPossible (where there some that didn't match?!?), pT, eta, phi, X0, Y0, Z0,
//    also match to the silicone tracklets

class PHCompositeNode;
/* class TNtuple; */


class PHG4TruthInfoContainer;
class SvtxTrackMap;
class TrackSeedContainer;
class TrkrClusterContainer;
class TrkrTruthTrackContainer;

// truth information to match to the track
/* class TruthRecoTrackMatching : public Fun4AllBase { */

class TruthRecoTrackMatching : public SubsysReco 
{
  public:
  /* TruthRecoTrackMatching(const std::string &name = "TruthRecoTrackMatching"); */
  TruthRecoTrackMatching() {};
  ~TruthRecoTrackMatching() override = default;
  int Init(PHCompositeNode          *) override { return 0; };
  int InitRun(PHCompositeNode       *) override;
  int process_event(PHCompositeNode *) override;
  int End(PHCompositeNode           *) override;


  private:
  int createNodes(PHCompositeNode *topNode);

  // INPUT
  PHG4TruthInfoContainer  *m_PHG4TruthInfoContainer  {nullptr}; // Get the truth track ids
  SvtxTrackMap            *m_SvtxTrackMap            {nullptr}; 
  TrkrClusterContainer    *m_TruthClusterContainer   {nullptr};
  TrkrClusterContainer    *m_RecoClusterContainer    {nullptr};
  TrkrTruthTrackContainer *m_TrkrTruthTrackContainer {nullptr};

  // The sorting can get:
  // 0. EmbedRecoPairs
  // 1. EmbedNotMatched
  // 2. RecoNotMatched
  
  // locally - maybe?
  TrackSeedContainer      *m_tpcTrackSeedContainer   {nullptr}; // Get the seeds from the tracks to get the clusters
  // OUTPUT

  /*
   * Files:
   ~/coresoftware/simulation/g4simulation/g4main/PHG4TruthInfoContainer.h
   ~/coresoftware/offline/packages/trackbase_historic/SvtxTrackMap_v2.h
   ~/coresoftware/offline/packages/trackreco/PHTrackCleaner.h
   */

  // The low down:
  //
  // Input:
  // SvtxTrackMap_v2 ``name?'' -> Container for SvtxTrack* (I will pull SvtxTrack_v4 from them)
  //   ┕▶ SvtexTrack_v4   -> Get from SvtxTrackMap; can loop through them to find matched
  // TrackSeedContainer_v1 ``name?''
  //
  //
  // Output:
  // m_TrkrTruthTrackContainer ->
  //
  //
  /*
  // ----------------------
  // ----------------------
  // ----------------------
  // ----------------------
  // ----------------------
  // ----------------------
  // Match to:
  //    Tracks - pT, phi, eta -- > Where from? Are these the track seeds? Which node, which name
  from trackbase_historic --  SvtxTrack_v4(); -> Has pointer to the Silicon and TPC seed
    The seeds have the cluster keys, pt, eta, and phi for the matching
      for best track, look to the track seeds and loop over
    These tracks are already cleaned by track cleaner (after ACTS fit)

  SvtxTrackSeed_v1.h  -> Has the track seeds' pointers

  For example see the PH..TrackFitter, and in the trackreco/PHTrackCleaner.cc
    (shows get nodes etc...) (for now just the TPC seed cluster, but ultimately will want
        the silicon clusters, too)

  All the truth tracks are PHG4Particle* objects, on _______ node from PHG4TruthInfoContainer::Map 
  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo") :: ParticleMap
  // so loop through  PHG4Truth

  
  // ----------------------
  // ----------------------
  // ----------------------
  // ----------------------
  // ----------------------
  */
  //---------------------- SvtxTrack

  // SubsysReco provides:
  //   /// Called at the end of all processing.
  //   virtual int End(PHCompositeNode * /*topNode*/) { return 0; }

  //   /// Called at the end of each run.
  //   virtual int EndRun(const int /*runnumber*/) { return 0; }

  //   /** Called during initialization.
  //       Typically this is where you can book histograms, and e.g.
  //       register them to Fun4AllServer (so they can be output to file
  //       using Fun4AllServer::dumpHistos() method).
  //    */
  //   virtual int Init(PHCompositeNode * /*topNode*/) { return 0; }

  //   /** Called for first event when run number is known.
  //       Typically this is where you may want to fetch data from
  //       database, because you know the run number.
  //    */
  //   virtual int InitRun(PHCompositeNode * /*topNode*/) { return 0; }

  //   /** Called for each event.
  //       This is where you do the real work.
  //   */
  //   virtual int process_event(PHCompositeNode * /*topNode*/) { return 0; }

  //   /// Reset.
  //   virtual int Reset(PHCompositeNode * /*topNode*/) { return 0; }

  //   /// Clean up after each event.
  //   virtual int ResetEvent(PHCompositeNode * /*topNode*/) { return 0; }

  //   void Print(const std::string & /*what*/ = "ALL") const override {}
};


#endif
