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
class ActsPropagator;

using Trajectory = ActsExamples::Trajectories;

class MakeMilleFiles : public SubsysReco
{
 public:
  MakeMilleFiles(const std::string& name = "MakeMilleFiles");

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  void set_binary(bool bin) { _binary = bin; }
  void set_constraintfile_name(const std::string& file) { m_constraintFileName = file; }
  void set_datafile_name(const std::string& file) { data_outfilename = file; }
  void set_steeringfile_name(const std::string& file) { steering_outfilename = file; }
  void set_mvtx_grouping(int group) { mvtx_group = (AlignmentDefs::mvtxGrp) group; }
  void set_intt_grouping(int group) { intt_group = (AlignmentDefs::inttGrp) group; }
  void set_tpc_grouping(int group) { tpc_group = (AlignmentDefs::tpcGrp) group; }
  void set_mms_grouping(int group) { mms_group = (AlignmentDefs::mmsGrp) group; }
  void use_event_vertex(bool useit) { m_useEventVertex = useit; }
  void set_layer_fixed(unsigned int layer);
  void set_mvtx_layer_fixed(unsigned int layer, unsigned int clamshell);
  void set_layer_gparam_fixed(unsigned int layer, unsigned int param);
  void set_layer_lparam_fixed(unsigned int layer, unsigned int param);

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
  void set_vtx_sigma(float xysig, float zsig)
  {
    m_vtxSigma(0) = xysig;
    m_vtxSigma(1) = zsig;
  }

 private:
  Mille* _mille;

  int GetNodes(PHCompositeNode* topNode);
  Acts::Vector3 getEventVertex();

  bool is_layer_fixed(unsigned int layer);

  bool is_layer_param_fixed(unsigned int layer, unsigned int param, std::set<std::pair<unsigned int, unsigned int>>& param_fixed);

  bool is_tpc_sector_fixed(unsigned int layer, unsigned int sector, unsigned int side);
  bool is_mvtx_layer_fixed(unsigned int layer, unsigned int clamshell);
  void addTrackToMilleFile(SvtxAlignmentStateMap::StateVec& statevec);
  void getGlobalVtxDerivativesXY(SvtxTrack* track,
                                 const Acts::Vector3& vertex,
                                 float glblvtx_derivative[SvtxAlignmentState::NRES][3]);
  bool getLocalVtxDerivativesXY(SvtxTrack* track,
                                ActsPropagator& propagator,
                                const Acts::Vector3& vertex,
                                float lclvtx_derivative[SvtxAlignmentState::NRES][SvtxAlignmentState::NLOC]);
  Acts::Vector3 localToGlobalVertex(SvtxTrack* track,
                                    const Acts::Vector3& vertex,
                                    const Acts::Vector3& localx) const;
  void getProjectionVtxXY(SvtxTrack* track, const Acts::Vector3& vertex,
                          Acts::Vector3& projx, Acts::Vector3& projy);

  std::string data_outfilename = ("mille_output_data_file.bin");
  std::string steering_outfilename = ("steer.txt");

  bool m_useEventVertex = false;
  bool _binary = true;

  Acts::Vector2 m_vtxSigma = {0.1, 0.1};

  std::map<unsigned int, float> m_layerMisalignment;
  std::set<unsigned int> fixed_sectors;
  // set default groups to lowest level
  AlignmentDefs::mvtxGrp mvtx_group = AlignmentDefs::mvtxGrp::snsr;
  AlignmentDefs::inttGrp intt_group = AlignmentDefs::inttGrp::chp;
  AlignmentDefs::tpcGrp tpc_group = AlignmentDefs::tpcGrp::htst;
  AlignmentDefs::mmsGrp mms_group = AlignmentDefs::mmsGrp::tl;

  std::set<unsigned int> fixed_layers;
  std::set<std::pair<unsigned int, unsigned int>> fixed_mvtx_layers;
  std::set<std::pair<unsigned int, unsigned int>> fixed_layer_gparams, fixed_layer_lparams;

  std::string m_constraintFileName = "mp2con.txt";
  std::ofstream m_constraintFile;

  SvtxTrackMap* _track_map{nullptr};
  SvtxAlignmentStateMap* _state_map{nullptr};
  ActsGeometry* _tGeometry{nullptr};
  TrkrClusterContainer* _cluster_map{nullptr};
  ClusterErrorPara _ClusErrPara;
};

#endif  // MAKEMILLEFILES_H
