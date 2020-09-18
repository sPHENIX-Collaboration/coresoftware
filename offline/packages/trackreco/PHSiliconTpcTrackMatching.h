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
class TF1;

#include <trackreco/PHTrackPropagating.h>

class PHSiliconTpcTrackMatching : public PHTrackPropagating
{
 public:

  PHSiliconTpcTrackMatching(const std::string &name = "PHSiliconTpcTrackMatching");

  virtual ~PHSiliconTpcTrackMatching();

  void set_track_map_name_silicon(const std::string &map_name) { _track_map_name_silicon = map_name; }
  void set_phi_search_window(const double win){_phi_search_win = win;}
  void set_eta_search_window(const double win){_eta_search_win = win;}
  void set_search_par_values(const double p0, const double p1, const double p2){   ;}

 protected:
  int Setup(PHCompositeNode* topNode) override;

  int Process() override;

  int End() override;
  
 private:

  int GetNodes(PHCompositeNode* topNode);

  std::string _track_map_name_silicon;
  double _phi_search_win = 0.03;
  double _eta_search_win = 0.02;
  
  SvtxTrackMap *_track_map_silicon{nullptr};
  SvtxTrack *_tracklet_tpc{nullptr};
  SvtxTrack *_tracklet_si{nullptr};

  TF1 *fdphi{nullptr};
  double par0 =  -0.000650;
  double par1 =  0.13373;
  double par2 =  0.98298;

};

#endif // PHTRUTHSILICONASSOCIATION_H
