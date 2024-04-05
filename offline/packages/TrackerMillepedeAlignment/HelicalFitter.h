// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef HELICALFITTER_H
#define HELICALFITTER_H

#include "AlignmentDefs.h"

#include <tpc/TpcClusterZCrossingCorrection.h>
#include <tpc/TpcDistortionCorrection.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrackFitUtils.h>

#include <phparameter/PHParameterInterface.h>

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>

class PHCompositeNode;
class TrackSeedContainer;
class TrackSeed;
class TrkrClusterContainer;
class TF1;
class TNtuple;
class TFile;
class TpcDistortionCorrectionContainer;
class Mille;
class SvtxTrackSeed;
class SvtxTrackMap;
class SvtxAlignmentStateMap;
class SvtxTrack;

class HelicalFitter : public SubsysReco, public PHParameterInterface
{
 public:
  HelicalFitter(const std::string& name = "HelicalFitter");

  ~HelicalFitter() override;

  void SetDefaultParameters() override;

  void set_field_dir(const double rescale)
  {
    _fieldDir = -1;
    if (rescale > 0)
    {
      _fieldDir = 1;
    }
  }
  void set_field(const std::string& field) { _field = field; }

  int InitRun(PHCompositeNode* topNode) override;

  int process_event(PHCompositeNode*) override;

  int End(PHCompositeNode*) override;

  void set_silicon_track_map_name(const std::string& map_name) { _silicon_track_map_name = map_name; }
  void set_track_map_name(const std::string& map_name) { _track_map_name = map_name; }

  void set_use_event_vertex(bool flag) { use_event_vertex = flag; }
  void set_datafile_name(const std::string& file) { data_outfilename = file; }
  void set_steeringfile_name(const std::string& file) { steering_outfilename = file; }
  void set_mvtx_grouping(int group) { mvtx_grp = (AlignmentDefs::mvtxGrp) group; }
  void set_intt_grouping(int group) { intt_grp = (AlignmentDefs::inttGrp) group; }
  void set_tpc_grouping(int group) { tpc_grp = (AlignmentDefs::tpcGrp) group; }
  void set_mms_grouping(int group) { mms_grp = (AlignmentDefs::mmsGrp) group; }
  void set_test_output(bool test) { test_output = test; }
  void set_intt_layer_fixed(unsigned int layer);
  void set_mvtx_layer_fixed(unsigned int layer, unsigned int clamshell);
  void set_tpc_sector_fixed(unsigned int region, unsigned int sector, unsigned int side);
  void set_layer_param_fixed(unsigned int layer, unsigned int param);
  void set_ntuplefile_name(const std::string& file) { ntuple_outfilename = file; }

  void set_fitted_subsystems(bool si, bool tpc, bool full)
  {
    fitsilicon = si;
    fittpc = tpc;
    fitfulltrack = full;
  }

  void set_error_inflation_factor(unsigned int layer, float factor)
  {
    _layerMisalignment.insert(std::make_pair(layer, factor));
  }

  void set_vertex(Acts::Vector3 inVertex, Acts::Vector3 inVertex_uncertainty)
  {
    vertexPosition = inVertex;
    vertexPosUncertainty = inVertex_uncertainty;
  }

  void set_vtx_sigma(float x_sigma, float y_sigma)
  {
    vtx_sigma(0) = x_sigma;
    vtx_sigma(1) = y_sigma;
  }

  // utility functions for analysis modules
  std::vector<float> fitClusters(std::vector<Acts::Vector3>& global_vec, std::vector<TrkrDefs::cluskey> cluskey_vec);
  void getTrackletClusters(TrackSeed* _track, std::vector<Acts::Vector3>& global_vec, std::vector<TrkrDefs::cluskey>& cluskey_vec);
  Acts::Vector3 get_helix_pca(std::vector<float>& fitpars, const Acts::Vector3& global);
  void correctTpcGlobalPositions(std::vector<Acts::Vector3> global_vec, std::vector<TrkrDefs::cluskey> cluskey_vec);
  unsigned int addSiliconClusters(std::vector<float>& fitpars, std::vector<Acts::Vector3>& global_vec, std::vector<TrkrDefs::cluskey>& cluskey_vec);

  void set_dca_cut(float dca) { dca_cut = dca; }

 private:
  Mille* _mille;

  int GetNodes(PHCompositeNode* topNode);
  int CreateNodes(PHCompositeNode* topNode);
  void getTrackletClusterList(TrackSeed* tracklet, std::vector<TrkrDefs::cluskey>& cluskey_vec);

  Acts::Vector3 getPCALinePoint(const Acts::Vector3& global, const Acts::Vector3& tangent, const Acts::Vector3& posref);
  Acts::Vector3 get_line_plane_intersection(const Acts::Vector3& PCA, const Acts::Vector3& tangent,
                                            const Acts::Vector3& sensor_center, const Acts::Vector3& sensor_normal);
  std::pair<Acts::Vector3, Acts::Vector3> get_helix_tangent(const std::vector<float>& fitpars, Acts::Vector3 global);
  Acts::Vector3 get_helix_surface_intersection(const Surface& surf, std::vector<float>& fitpars, Acts::Vector3 global);
  Acts::Vector3 get_helix_surface_intersection(const Surface& surf, std::vector<float>& fitpars, Acts::Vector3 global, Acts::Vector3& pca, Acts::Vector3& tangent);

  Acts::Vector3 get_helix_vtx(Acts::Vector3 event_vtx, const std::vector<float>& fitpars);

  float convertTimeToZ(TrkrDefs::cluskey cluster_key, TrkrCluster* cluster);
  void makeTpcGlobalCorrections(TrkrDefs::cluskey cluster_key, short int crossing, Acts::Vector3& global);

  Acts::Vector2 getClusterError(TrkrCluster* cluster, TrkrDefs::cluskey cluskey, Acts::Vector3& global);

  bool is_tpc_sector_fixed(unsigned int layer, unsigned int sector, unsigned int side);
  bool is_mvtx_layer_fixed(unsigned int layer, unsigned int stave);
  bool is_intt_layer_fixed(unsigned int layer);
  bool is_layer_param_fixed(unsigned int layer, unsigned int param);

  void getLocalDerivativesXY(const Surface& surf, const Acts::Vector3& global, const std::vector<float>& fitpars, float lcl_derivativeX[5], float lcl_derivativeY[5], unsigned int layer);

  void getLocalVtxDerivativesXY(SvtxTrack& track, const Acts::Vector3& track_vtx, const std::vector<float>& fitpars, float lcl_derivativeX[5], float lcl_derivativeY[5]);

  void getGlobalDerivativesXY(const Surface& surf, Acts::Vector3 global, const Acts::Vector3& fitpoint, const std::vector<float>& fitpars, float glb_derivativeX[6], float glbl_derivativeY[6], unsigned int layer);

  void getGlobalVtxDerivativesXY(SvtxTrack& track, const Acts::Vector3& track_vtx, float glbl_derivativeX[3], float glbl_derivativeY[3]);

  void get_projectionXY(const Surface& surf, const std::pair<Acts::Vector3, Acts::Vector3>& tangent, Acts::Vector3& projX, Acts::Vector3& projY);
  void get_projectionVtxXY(SvtxTrack& track, const Acts::Vector3& event_vtx, Acts::Vector3& projX, Acts::Vector3& projY);

  float getVertexResidual(Acts::Vector3 vtx);

  void get_dca(SvtxTrack& track, float& dca3dxy, float& dca3dz, float& dca3dxysigma, float& dca3dzsigma, const Acts::Vector3& vertex);
  Acts::Vector3 globalvtxToLocalvtx(SvtxTrack& track, const Acts::Vector3& event_vertex);
  Acts::Vector3 globalvtxToLocalvtx(SvtxTrack& track, const Acts::Vector3& event_vertex, Acts::Vector3 PCA);
  Acts::Vector3 localvtxToGlobalvtx(SvtxTrack& track, const Acts::Vector3& event_vtx, const Acts::Vector3& PCA);

  TpcClusterZCrossingCorrection m_clusterCrossingCorrection;
  TpcDistortionCorrectionContainer* _dcc_static{nullptr};
  TpcDistortionCorrectionContainer* _dcc_average{nullptr};
  TpcDistortionCorrectionContainer* _dcc_fluctuation{nullptr};

  bool test_output = false;

  ClusterErrorPara _ClusErrPara;

  std::set<std::pair<unsigned int, unsigned int>> fixed_mvtx_layers;
  std::set<unsigned int> fixed_intt_layers;
  std::set<unsigned int> fixed_sectors;
  std::set<std::pair<unsigned int, unsigned int>> fixed_layer_params;

  // set default groups to lowest level
  AlignmentDefs::mvtxGrp mvtx_grp = AlignmentDefs::mvtxGrp::snsr;
  AlignmentDefs::inttGrp intt_grp = AlignmentDefs::inttGrp::chp;
  AlignmentDefs::tpcGrp tpc_grp = AlignmentDefs::tpcGrp::htst;
  AlignmentDefs::mmsGrp mms_grp = AlignmentDefs::mmsGrp::tl;

  /// tpc distortion correction utility class
  TpcDistortionCorrection _distortionCorrection;

  //  TrackSeedContainer *_svtx_seed_map{nullptr};
  TrackSeedContainer* _track_map_tpc{nullptr};
  TrackSeedContainer* _track_map_silicon{nullptr};
  TrkrClusterContainer* _cluster_map{nullptr};
  ActsGeometry* _tGeometry{nullptr};

  std::string data_outfilename{"mille_helical_output_data_file.bin"};
  std::string steering_outfilename{"steer_helical.txt"};
  std::string ntuple_outfilename{"HF_ntuple.root"};

  bool fitsilicon{true};
  bool fittpc{false};
  bool fitfulltrack{false};

  float dca_cut{0.19};  // cm

  SvtxTrackMap* m_trackmap{nullptr};
  SvtxAlignmentStateMap* m_alignmentmap{nullptr};

  std::string _field;
  int _fieldDir{-1};
  std::map<unsigned int, float> _layerMisalignment;

  std::string _track_map_name{"TpcTrackSeedContainer"};
  std::string _silicon_track_map_name{"SiliconTrackSeedContainer"};

  bool make_ntuple{true};
  TNtuple* ntp{nullptr};
  TNtuple* track_ntp{nullptr};
  TFile* fout{nullptr};

  bool use_event_vertex{false};

  int event{0};

  Acts::Vector3 vertexPosition;
  Acts::Vector3 vertexPosUncertainty;
  Acts::Vector2 vtx_sigma;
};

#endif  // HELICALFITTER_H
