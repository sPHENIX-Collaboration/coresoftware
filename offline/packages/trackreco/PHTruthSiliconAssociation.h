// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHTRUTHSILICONASSOCIATION_H
#define PHTRUTHSILICONASSOCIATION_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>
#include <set>

class PHCompositeNode;
class SvtxTrackMap;
class SvtxTrack;
class SvtxVertexMap;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;
class PHG4TruthInfoContainer;
class PHG4HitContainer;
class PHG4Particle;
class AssocInfoContainer;


class PHTruthSiliconAssociation : public SubsysReco
{
 public:

  PHTruthSiliconAssociation(const std::string &name = "PHTruthSiliconAssociation");

  virtual ~PHTruthSiliconAssociation();

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;
  
 private:

  int GetNodes(PHCompositeNode* topNode);

  PHG4Particle* getG4PrimaryParticle(SvtxTrack *track);
  std::set<TrkrDefs::cluskey> getSiliconClustersFromParticle(PHG4Particle* g4particle);
  
  PHG4TruthInfoContainer *_truthInfo{nullptr};
  PHG4HitContainer *_g4hits_tpc{nullptr};
  PHG4HitContainer *_g4hits_mvtx{nullptr};
  PHG4HitContainer *_g4hits_intt{nullptr};
  
  TrkrClusterContainer *_cluster_map{nullptr};
  TrkrClusterHitAssoc *_cluster_hit_map{nullptr};
  TrkrHitTruthAssoc *_hit_truth_map{nullptr};
  SvtxTrackMap *_track_map{nullptr};
  AssocInfoContainer *_assoc_container{nullptr};
  SvtxTrack *_tracklet{nullptr};
  SvtxVertexMap * _vertex_map{nullptr};
  PHG4TruthInfoContainer *_truthinfo{nullptr};

};

#endif // PHTRUTHSILICONASSOCIATION_H
