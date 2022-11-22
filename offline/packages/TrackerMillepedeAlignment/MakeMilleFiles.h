// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 *  \file		  MakeMilleFiles.h
 *  \brief		Class for moving corrected TPC clusters to the nearest TPC readout layer radius
 *  \author	 Tony Frawley <afrawley@fsu.edu>
 */

#ifndef MAKEMILLEFILES_H
#define MAKEMILLEFILES_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/ActsTransformations.h>

#include <trackbase_historic/SvtxTrackState_v1.h>

#include <tpc/TpcClusterZCrossingCorrection.h>
#include <tpc/TpcDistortionCorrection.h>
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <ActsExamples/EventData/Track.hpp>
#include <ActsExamples/EventData/Trajectories.hpp>

class PHCompositeNode;
class PHG4TpcCylinderGeomContainer;
class SvtxTrack;
class SvtxTrackMap;
class TrkrCluster;
class TrkrClusterContainer;
class TpcDistortionCorrectionContainer;
class ClusterErrorPara;
class Mille;

using Trajectory = ActsExamples::Trajectories;

enum siliconGroup
{
  sensor,
  stave,
  barrel
};
enum tpcGroup
{
  hitset,
  sector,
  tpc
};
enum mmsGroup
{
  tile,
  mms
};

/*
 * Class which contains alignment information to be written to pede, 
 * obtained from Acts::Trajectories and track states
 * Also defines all matrix and residual dimensions
 */
class AlignmentState
{
 public:
  /// The number of global (alignment) parameters
  static const int NGL = 6;
  /// The number of local (track state) parameters
  static const int NLC = 6;
  /// The number of residuals per state (e.g. 2D or 3D)
  static const int NRES = 3;

  using GlobalMatrix = Acts::ActsMatrix<NRES, NGL>;
  using LocalMatrix = Acts::ActsMatrix<NRES, NLC>;
  using ResidualVector = Eigen::Matrix<Acts::ActsScalar, NRES, 1>;

  AlignmentState(size_t index, ResidualVector res,
                 GlobalMatrix ralign,
                 LocalMatrix rtrack,
                 Acts::Vector3 clusglob)
    : m_tsIndex(index)
    , m_residual(res)
    , m_dResAlignmentPar(ralign)
    , m_dResTrackPar(rtrack)
    , m_clusglob(clusglob)
  {
  }

  void set_residual(const ResidualVector& res) { m_residual = res; }
  void set_dResAlignmentPar(const GlobalMatrix& d)
  {
    m_dResAlignmentPar = d;
  }
  void set_dResTrackPar(const LocalMatrix& d)
  {
    m_dResTrackPar = d;
  }

  const ResidualVector& get_residual() const { return m_residual; }
  const GlobalMatrix& get_dResAlignmentPar() const
  {
    return m_dResAlignmentPar;
  }
  const LocalMatrix& get_dResTrackPar() const
  {
    return m_dResTrackPar;
  }
  void set_tsIndex(const size_t index) { m_tsIndex = index; }
  const size_t& get_tsIndex() const { return m_tsIndex; }
  const Acts::Vector3& get_clusglob() const { return m_clusglob; }

 private:
  size_t m_tsIndex;
  ResidualVector m_residual;
  GlobalMatrix m_dResAlignmentPar;
  LocalMatrix m_dResTrackPar;
  Acts::Vector3 m_clusglob;
};

using AlignmentStateMap = std::map<TrkrDefs::cluskey, AlignmentState>;

class MakeMilleFiles : public SubsysReco
{
 public:
  MakeMilleFiles(const std::string& name = "MakeMilleFiles");

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  void set_binary(bool bin) { _binary = bin; }

  void set_datafile_name(const std::string& file) { data_outfilename = file; }
  void set_steeringfile_name(const std::string& file) { steering_outfilename = file; }
  void set_silicon_grouping(int group) { si_group = (siliconGroup) group; }
  void set_tpc_grouping(int group) { tpc_group = (tpcGroup) group; }
  void set_mms_grouping(int group) { mms_group = (mmsGroup) group; }

 private:
  Mille* _mille;

  std::map<const unsigned int, Trajectory>* _trajectories;

  int GetNodes(PHCompositeNode* topNode);
  Acts::Vector3 getPCALinePoint(Acts::Vector3 global, SvtxTrackState* state);
  std::vector<Acts::Vector3> getDerivativesAlignmentAngles(Acts::Vector3& global,
                                                           TrkrDefs::cluskey cluster_key,
                                                           TrkrCluster* cluster,
                                                           Surface surface, int crossing);
  SvtxTrack::StateIter getStateIter(Acts::Vector3& global, SvtxTrack* track);
  void makeTpcGlobalCorrections(TrkrDefs::cluskey cluster_key,
                                short int crossing, Acts::Vector3& global);
  float convertTimeToZ(TrkrDefs::cluskey cluster_key, TrkrCluster* cluster);
  Acts::Transform3 makePerturbationTransformation(Acts::Vector3 angles);
  int getLabelBase(Acts::GeometryIdentifier id);
  int getTpcRegion(int layer);

  AlignmentStateMap getAlignmentStates(const Trajectory& traj,
                                       SvtxTrack* track, short int crossing);
  void addTrackToMilleFile(AlignmentStateMap& alignStates, const Trajectory& traj);

  std::map<int, float> derivativeGL;
  std::string data_outfilename = ("mille_output_data_file.bin");
  std::string steering_outfilename = ("steer.txt");

  /// tpc distortion correction utility class
  TpcDistortionCorrection _distortionCorrection;
  bool _binary = true;
  unsigned int _cluster_version = 4;

  bool m_useAnalytic = true;

  ClusterErrorPara _ClusErrPara;

  float sensorAngles[3] = {0.1, 0.1, 0.2};  // perturbation values for each alignment angle

  // set default groups to lowest level
  siliconGroup si_group = siliconGroup::sensor;
  tpcGroup tpc_group = tpcGroup::hitset;
  mmsGroup mms_group = mmsGroup::tile;

  int nsensors_stave[7] = {9,9,9,4,4,4,4};

  std::map<unsigned int, unsigned int> base_layer_map = {{10, 0}, {12, 3}, {14, 7}, {16, 55}};

  SvtxTrackMap* _track_map{nullptr};
  TrkrClusterContainer* _cluster_map{nullptr};
  ActsGeometry* _tGeometry{nullptr};

  TpcClusterZCrossingCorrection m_clusterCrossingCorrection;
  TpcDistortionCorrectionContainer* _dcc_static{nullptr};
  TpcDistortionCorrectionContainer* _dcc_average{nullptr};
  TpcDistortionCorrectionContainer* _dcc_fluctuation{nullptr};
};

#endif  // MAKEMILLEFILES_H
