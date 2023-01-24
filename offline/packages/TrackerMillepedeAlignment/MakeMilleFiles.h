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
#include <trackbase_historic/SvtxAlignmentStateMap.h>

#include <ActsExamples/EventData/Trajectories.hpp>
class PHCompositeNode;
class PHG4TpcCylinderGeomContainer;
class SvtxTrack;
class SvtxTrackMap;
class TrkrCluster;
class TrkrClusterContainer;
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
  void set_layer_fixed(unsigned int layer);
  void set_layer_param_fixed(unsigned int layer, unsigned int param);
  void set_cluster_version(unsigned int v) { _cluster_version = v; }

 private:
  Mille* _mille;

  int GetNodes(PHCompositeNode* topNode);
  Acts::Vector3 getPCALinePoint(Acts::Vector3 global, SvtxTrackState* state);
  std::vector<Acts::Vector3> getDerivativesAlignmentAngles(Acts::Vector3& global,
                                                           TrkrDefs::cluskey cluster_key,
                                                           TrkrCluster* cluster,
                                                           Surface surface, int crossing);

  int getLabelBase(Acts::GeometryIdentifier id);
  int getTpcRegion(int layer);

  bool is_layer_fixed(unsigned int layer);
  bool is_layer_param_fixed(unsigned int layer, unsigned int param);
  void printBuffers(int index, Acts::Vector3 residual, Acts::Vector3 clus_sigma, float lcl_derivative[], float glbl_derivative[], int glbl_label[]);

  void addTrackToMilleFile(SvtxAlignmentStateMap::StateVec statevec);

  std::map<int, float> derivativeGL;
  std::string data_outfilename = ("mille_output_data_file.bin");
  std::string steering_outfilename = ("steer.txt");

  bool _binary = true;
  unsigned int _cluster_version = 4;

  bool m_useAnalytic = true;

  // set default groups to lowest level
  siliconGroup si_group = siliconGroup::sensor;
  tpcGroup tpc_group = tpcGroup::hitset;
  mmsGroup mms_group = mmsGroup::tile;

  int nsensors_stave[7] = {9, 9, 9, 4, 4, 4, 4};

  std::set<unsigned int> fixed_layers;
  std::set<std::pair<unsigned int, unsigned int>> fixed_layer_params;

  std::map<unsigned int, unsigned int> base_layer_map = {{10, 0}, {12, 3}, {14, 7}, {16, 55}};

  SvtxTrackMap* _track_map{nullptr};
  SvtxAlignmentStateMap* _state_map{nullptr};
  ActsGeometry* _tGeometry{nullptr};
  TrkrClusterContainer* _cluster_map{nullptr};
  ClusterErrorPara _ClusErrPara;
};

#endif  // MAKEMILLEFILES_H
