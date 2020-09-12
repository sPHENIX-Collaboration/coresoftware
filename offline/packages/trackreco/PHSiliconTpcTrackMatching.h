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
  
 private:

  int GetNodes(PHCompositeNode* topNode);

  std::string _track_map_name_silicon;
  double _phi_search_win = 0.05;
  double _eta_search_win = 0.02;
  
  SvtxTrackMap *_track_map_silicon{nullptr};
  SvtxTrack *_tracklet_tpc{nullptr};
  SvtxTrack *_tracklet_si{nullptr};
};

#endif // PHTRUTHSILICONASSOCIATION_H
