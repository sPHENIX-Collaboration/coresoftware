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
class SvtxTrackSeed;

enum siliconGrp {snsr, stv, brrl};
enum tpcGrp {htst, sctr, tp};
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
  void set_layer_fixed(unsigned int layer);
  void set_layer_param_fixed(unsigned int layer, unsigned int param);
  void set_cluster_version(unsigned int v) { _cluster_version = v; }

 private:

Mille* _mille;

  int GetNodes(PHCompositeNode* topNode);

  Acts::Vector3 get_helix_pca(std::vector<float>& fitpars, Acts::Vector3 global);
  Acts::Vector3 getPCALinePoint(Acts::Vector3 global, Acts::Vector3 tangent, Acts::Vector3 posref);
  Acts::Vector2 get_circle_point_pca(float radius, float x0, float y0, Acts::Vector3 global);

  int getLabelBase(Acts::GeometryIdentifier id);
  Acts::Transform3 makePerturbationTransformation(Acts::Vector3 angles);
  std::vector<Acts::Vector3> getDerivativesAlignmentAngles(Acts::Vector3& global, TrkrDefs::cluskey cluster_key, TrkrCluster* cluster, Surface surface, int crossing);
  float convertTimeToZ(TrkrDefs::cluskey cluster_key, TrkrCluster *cluster);
  void makeTpcGlobalCorrections(TrkrDefs::cluskey cluster_key, short int crossing, Acts::Vector3& global);
  int getTpcRegion(int layer);

  void getTrackletClusters(TrackSeed *_track, std::vector<Acts::Vector3>& global_vec, std::vector<TrkrDefs::cluskey>& cluskey_vec);
  std::vector<float> fitClusters(std::vector<Acts::Vector3>& global_vec, std::vector<TrkrDefs::cluskey> cluskey_vec);
  Acts::Vector3 getClusterError(TrkrCluster *cluster, TrkrDefs::cluskey cluskey, Acts::Vector3& global);
  void getGlobalLabels(Surface surf, int glbl_label[]);
  void getLocalDerivativesX(Acts::Vector3& pca, std::vector<float>& fitpars, float lcl_derivative[]);
  void getLocalDerivativesY(Acts::Vector3& pca, std::vector<float>& fitpars, float lcl_derivative[]);
  void getLocalDerivativesZ(Acts::Vector3& global, float lcl_derivative[]);
  void getGlobalDerivativesX( std::vector<Acts::Vector3> angleDerivs, float glbl_derivatives[], unsigned int layer);
  void getGlobalDerivativesY( std::vector<Acts::Vector3> angleDerivs, float glbl_derivatives[], unsigned int layer);
  void getGlobalDerivativesZ( std::vector<Acts::Vector3> angleDerivs, float glbl_derivatives[], unsigned int layer);
  void printBuffers(int index, Acts::Vector3 residual, Acts::Vector3 clus_sigma, float lcl_derivative[], float glbl_derivative[], int glbl_label[]);
  bool is_layer_fixed(unsigned int layer);
  bool is_layer_param_fixed(unsigned int layer, unsigned int param);
  unsigned int addSiliconClusters(std::vector<float>& fitpars, std::vector<Acts::Vector3>& global_vec,  std::vector<TrkrDefs::cluskey>& cluskey_vec);

  TpcClusterZCrossingCorrection m_clusterCrossingCorrection;
  TpcDistortionCorrectionContainer* _dcc_static{nullptr};
  TpcDistortionCorrectionContainer* _dcc_average{nullptr};
  TpcDistortionCorrectionContainer* _dcc_fluctuation{nullptr};

  unsigned int _cluster_version = 4;
  bool test_output = false;

  ClusterErrorPara _ClusErrPara;

  float sensorAngles[3] = {0.1, 0.1, 0.2};  // perturbation values for each alignment angle
 
  std::set<unsigned int> fixed_layers;
  std::set<std::pair<unsigned int,unsigned int>> fixed_layer_params;

  // set default groups to lowest level
  siliconGrp si_grp = siliconGrp::snsr;
  tpcGrp tpc_grp = tpcGrp::htst;
  mmsGrp mms_grp = mmsGrp::tl;

  int nsensors_stave[7] = {9,9,9,4,4,4,4};

  std::map<unsigned int, unsigned int> base_layer_map = { {10, 0}, {12,3}, {14,7}, {16,55} };

 /// tpc distortion correction utility class
  TpcDistortionCorrection _distortionCorrection;

  //  TrackSeedContainer *_svtx_seed_map{nullptr};
  TrackSeedContainer *_track_map_tpc{nullptr};
  TrackSeedContainer *_track_map_silicon{nullptr};
  TrkrClusterContainer *_cluster_map{nullptr};
  ActsGeometry *_tGeometry{nullptr};

  std::string  data_outfilename = ("mille_helical_output_data_file.bin");  
  std::string  steering_outfilename = ("steer_helical.txt");  

  static const int NLC = 5;
  static const int NGL = 6;

  bool fitsilicon = false;
  bool fittpc = true;
  bool fitfulltrack = true;

  float dca_cut = 0.1;  // 1 mm

  std::string _field;
  int _fieldDir = -1;

  std::string _track_map_name = "TpcTrackSeedContainer";
  std::string _silicon_track_map_name = "SiliconTrackSeedContainer";
};

#endif // HELICALFITTER_H
