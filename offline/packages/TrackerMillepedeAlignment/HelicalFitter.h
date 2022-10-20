// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef HELICALFITTER_H
#define HELICALFITTER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
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
class Mille;

enum siliconGrp {snsr, stv, brrl};
enum tpcGrp {subsrf, sctr, tp};
enum mmsGrp {tl, mm};

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

 void set_datafile_name(const std::string& file) { data_outfilename = file;}
  void set_steeringfile_name(const std::string& file) { steering_outfilename = file;}
  void set_silicon_grouping(int group) {si_grp = (siliconGrp) group;}
  void set_tpc_grouping(int group) {tpc_grp = (tpcGrp) group;}
  void set_mms_grouping(int group) {mms_grp = (mmsGrp) group;}
  void set_test_output(bool test) {test_output = test;}

 private:

Mille* _mille;

  int GetNodes(PHCompositeNode* topNode);

  Acts::Vector3 get_helix_pca(float radius, float zslope, Acts::Vector3 helix_center, Acts::Vector3 global);
  Acts::Vector3 getPCALinePoint(Acts::Vector3 global, Acts::Vector3 tangent, Acts::Vector3 posref);
  Acts::Vector2 get_circle_point_pca(float radius, Acts::Vector3 center, Acts::Vector3 global);

  int getLabelBase(Acts::GeometryIdentifier id);
  Acts::Transform3 makePerturbationTransformation(Acts::Vector3 angles);
  std::vector<Acts::Vector3> getDerivativesAlignmentAngles(Acts::Vector3& global, TrkrDefs::cluskey cluster_key, TrkrCluster* cluster, Surface surface, int crossing);
  float convertTimeToZ(TrkrDefs::cluskey cluster_key, TrkrCluster *cluster);
  void makeTpcGlobalCorrections(TrkrDefs::cluskey cluster_key, short int crossing, Acts::Vector3& global);

  TpcClusterZCrossingCorrection m_clusterCrossingCorrection;
  TpcDistortionCorrectionContainer* _dcc_static{nullptr};
  TpcDistortionCorrectionContainer* _dcc_average{nullptr};
  TpcDistortionCorrectionContainer* _dcc_fluctuation{nullptr};

  unsigned int _cluster_version = 3;
  bool test_output = false;

  ClusterErrorPara _ClusErrPara;

  float sensorAngles[3] = {0.1, 0.1, 0.2};  // perturbation values for each alignment angle

  // set default groups to lowest level
  siliconGrp si_grp = siliconGrp::snsr;
  tpcGrp tpc_grp = tpcGrp::subsrf;
  mmsGrp mms_grp = mmsGrp::tl;

  int nstaves[7] = {12,16,20,12,12,16,16};

  std::map<unsigned int, unsigned int> base_layer_map = { {10, 0}, {12,3}, {14,7}, {16,55} };

 /// tpc distortion correction utility class
  TpcDistortionCorrection _distortionCorrection;

  //  TrackSeedContainer *_svtx_seed_map{nullptr};
  TrackSeedContainer *_track_map_tpc{nullptr};
  TrackSeedContainer *_track_map_silicon{nullptr};
  TrackSeed *_tracklet_tpc{nullptr};
  TrackSeed *_tracklet_si{nullptr};
  TrkrClusterContainer *_cluster_map{nullptr};
  ActsGeometry *_tGeometry{nullptr};

  std::string  data_outfilename = ("mille_helical_output_data_file.bin");  
  std::string  steering_outfilename = ("steer_helical.txt");  

  bool fitsilicon = true;
  bool fittpc = false;

  std::string _field;
  int _fieldDir = -1;

  std::string _track_map_name = "TpcTrackSeedContainer";
  std::string _silicon_track_map_name = "SiliconTrackSeedContainer";
};

#endif // HELICALFITTER_H
