// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 *  \file		  MakeMilleFiles.h
 *  \brief		Class for moving corrected TPC clusters to the nearest TPC readout layer radius
 *  \author	 Tony Frawley <afrawley@fsu.edu>
 */

#ifndef MAKEMILLEFILES_H
#define MAKEMILLEFILES_H

#include "AlignmentDefs.h"

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
  void set_silicon_grouping(int group) { si_group = (AlignmentDefs::siliconGrp) group; }
  void set_tpc_grouping(int group) { tpc_group = (AlignmentDefs::tpcGrp) group; }
  void set_mms_grouping(int group) { mms_group = (AlignmentDefs::mmsGrp) group; }
  void set_layer_fixed(unsigned int layer);
  void set_layer_param_fixed(unsigned int layer, unsigned int param);
  void set_cluster_version(unsigned int v) { _cluster_version = v; }
  void set_layers_fixed(unsigned int minlayer, unsigned int maxlayer);
  void set_error_inflation_factor(unsigned int layer, float factor)
  {
    m_layerMisalignment.insert(std::make_pair(layer, factor));
  }
  void set_tpc_sector_fixed(unsigned int region, unsigned int sector, 
			    unsigned int side)
 {
   // make a combined subsector index
   unsigned int subsector = region * 24 + side * 12 + sector;
   fixed_sectors.insert(subsector);
 }

 private:
  Mille* _mille;

  int GetNodes(PHCompositeNode* topNode);
  Acts::Vector3 getPCALinePoint(Acts::Vector3 global, SvtxTrackState* state);
  std::vector<Acts::Vector3> getDerivativesAlignmentAngles(Acts::Vector3& global,
                                                           TrkrDefs::cluskey cluster_key,
                                                           TrkrCluster* cluster,
                                                           Surface surface, int crossing);

  bool is_layer_fixed(unsigned int layer);
  bool is_layer_param_fixed(unsigned int layer, unsigned int param);
  bool is_tpc_sector_fixed(unsigned int layer, unsigned int sector, unsigned int side);
  void addTrackToMilleFile(SvtxAlignmentStateMap::StateVec statevec);

  std::map<int, float> derivativeGL;
  std::string data_outfilename = ("mille_output_data_file.bin");
  std::string steering_outfilename = ("steer.txt");

  bool _binary = true;
  unsigned int _cluster_version = 5;

  std::map<unsigned int, float> m_layerMisalignment;
  std::set<unsigned int> fixed_sectors;
  // set default groups to lowest level
  AlignmentDefs::siliconGrp si_group = AlignmentDefs::siliconGrp::snsr;
  AlignmentDefs::tpcGrp tpc_group = AlignmentDefs::tpcGrp::htst;
  AlignmentDefs::mmsGrp mms_group = AlignmentDefs::mmsGrp::tl;

  std::set<unsigned int> fixed_layers;
  std::set<std::pair<unsigned int, unsigned int>> fixed_layer_params;

  SvtxTrackMap* _track_map{nullptr};
  SvtxAlignmentStateMap* _state_map{nullptr};
  ActsGeometry* _tGeometry{nullptr};
  TrkrClusterContainer* _cluster_map{nullptr};
  ClusterErrorPara _ClusErrPara;
};

#endif  // MAKEMILLEFILES_H
