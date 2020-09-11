// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHSILICONTPCTRACKMATCHING_H
#define PHSILICONTPCTRACKMATCHING_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>
#include <set>
#include <vector>

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

#include <trackreco/PHTrackPropagating.h>

class PHSiliconTpcTrackMatching : public PHTrackPropagating
{
 public:

  PHSiliconTpcTrackMatching(const std::string &name = "PHSiliconTpcTrackMatching");

  virtual ~PHSiliconTpcTrackMatching();

  void set_track_map_name_silicon(const std::string &map_name) { _track_map_name_silicon = map_name; }

 protected:
  int Setup(PHCompositeNode* topNode) override;

  int Process() override;

  int End() override;

  //int Init(PHCompositeNode *topNode) override;
  //int InitRun(PHCompositeNode *topNode) override;
  //int process_event(PHCompositeNode *topNode) override;
  //int ResetEvent(PHCompositeNode *topNode) override;
  //int EndRun(const int runnumber) override;
  //int End(PHCompositeNode *topNode) override;
  //int Reset(PHCompositeNode * /*topNode*/) override;
  //void Print(const std::string &what = "ALL") const override;
  
 private:

  int GetNodes(PHCompositeNode* topNode);

  //std::vector<PHG4Particle*> getG4PrimaryParticle(SvtxTrack *track);
  //std::set<TrkrDefs::cluskey> getSiliconClustersFromParticle(PHG4Particle* g4particle);
  
  //PHG4TruthInfoContainer *_truthinfo{nullptr};
  //PHG4HitContainer *_g4hits_tpc{nullptr};
  //PHG4HitContainer *_g4hits_mvtx{nullptr};
  //PHG4HitContainer *_g4hits_intt{nullptr};

  std::string _track_map_name_silicon;
  double _phi_search_win = 0.02;
  double _eta_search_win = 0.01;
  
  SvtxTrackMap *_track_map_silicon{nullptr};
  SvtxTrack *_tracklet_tpc{nullptr};
  SvtxTrack *_tracklet_si{nullptr};
};

#endif // PHTRUTHSILICONASSOCIATION_H
