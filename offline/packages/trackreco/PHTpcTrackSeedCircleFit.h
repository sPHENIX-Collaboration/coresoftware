// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHTPCTRACKSEEDCIRCLEFIT_H
#define PHTPCTRACKSEEDCIRCLEFIT_H

#include <fun4all/SubsysReco.h>
#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase/ActsTrackingGeometry.h>
#include <tpc/TpcDistortionCorrection.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TrackSeed;
class TrackSeedContainer;
class TrkrCluster;
class TF1;
class TrkrClusterContainer;

class PHTpcTrackSeedCircleFit : public SubsysReco
{
 public:

  PHTpcTrackSeedCircleFit(const std::string &name = "PHTpcTrackSeedCircleFit");

  ~PHTpcTrackSeedCircleFit() override = default;

  int InitRun(PHCompositeNode* topNode) override;

  int process_event(PHCompositeNode*) override;

  int End(PHCompositeNode*) override;

  void use_truth_clusters(bool truth)
  { _use_truth_clusters = truth; }

  void set_track_map_name(const std::string &map_name) { _track_map_name = map_name; }
  void SetIteration(int iter){_n_iteration = iter;} 

 private:

  int GetNodes(PHCompositeNode* topNode);

  Acts::Vector3 getGlobalPosition( TrkrDefs::cluskey, TrkrCluster* cluster ) const;
						    
  ActsSurfaceMaps *_surfmaps{nullptr};
  ActsTrackingGeometry *_tGeometry{nullptr};
  TrackSeedContainer *_track_map{nullptr};
  
  bool _use_truth_clusters = false;
  TrkrClusterContainer *_cluster_map = nullptr;
  /// distortion correction container
  TpcDistortionCorrectionContainer* _dcc = nullptr;
 /// tpc distortion correction utility class
  TpcDistortionCorrection _distortionCorrection;

  int _n_iteration = 0;
  std::string _track_map_name = "SvtxTrackMap";

};

#endif // PHTRACKSEEDVERTEXASSOCIATION_H
