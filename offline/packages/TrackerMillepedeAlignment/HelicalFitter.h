// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef HELICALFITTER_H
#define HELICALFITTER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/ActsGeometry.h>
#include <phparameter/PHParameterInterface.h>

#include <tpc/TpcClusterZCrossingCorrection.h>
#include <tpc/TpcDistortionCorrection.h>

#include <string>
#include <map>

class PHCompositeNode;
class TrackSeedContainer;
class TrackSeed;
class TrkrClusterContainer;
class TF1;
class TpcDistortionCorrectionContainer;


class HelicalFitter : public SubsysReco, public PHParameterInterface
{
 public:

  HelicalFitter(const std::string &name = "HelicalFitter");

  ~HelicalFitter() override;

 void SetDefaultParameters() override;

  void set_field_dir(const double rescale)
  {
    _fieldDir = -1;
    if(rescale > 0)
      _fieldDir = 1;     
  }
  void set_field(const std::string &field) { _field = field;}

  int InitRun(PHCompositeNode* topNode) override;

  int process_event(PHCompositeNode*) override;

  int End(PHCompositeNode*) override;

  void set_silicon_track_map_name(const std::string &map_name) { _silicon_track_map_name = map_name; }
  void set_track_map_name(const std::string &map_name) { _track_map_name = map_name; }

 private:

  int GetNodes(PHCompositeNode* topNode);

  Acts::Vector3 get_helix_pca(float radius, float zslope, Acts::Vector3 helix_center, Acts::Vector3 global);
  Acts::Vector3 getPCALinePoint(Acts::Vector3 global, Acts::Vector3 tangent, Acts::Vector3 posref);
  Acts::Vector2 get_circle_point_pca(float radius, Acts::Vector3 center, Acts::Vector3 global);

  TpcClusterZCrossingCorrection m_clusterCrossingCorrection;
  TpcDistortionCorrectionContainer* _dcc_static{nullptr};
  TpcDistortionCorrectionContainer* _dcc_average{nullptr};
  TpcDistortionCorrectionContainer* _dcc_fluctuation{nullptr};

 /// tpc distortion correction utility class
  TpcDistortionCorrection _distortionCorrection;

  //  TrackSeedContainer *_svtx_seed_map{nullptr};
  TrackSeedContainer *_track_map_tpc{nullptr};
  TrackSeedContainer *_track_map_silicon{nullptr};
  TrackSeed *_tracklet_tpc{nullptr};
  TrackSeed *_tracklet_si{nullptr};
  TrkrClusterContainer *_cluster_map{nullptr};
  ActsGeometry *_tGeometry{nullptr};

  bool fitsilicon = true;
  bool fittpc = false;

  std::string _field;
  int _fieldDir = -1;

  std::string _track_map_name = "TpcTrackSeedContainer";
  std::string _silicon_track_map_name = "SiliconTrackSeedContainer";
};

#endif // HELICALFITTER_H
