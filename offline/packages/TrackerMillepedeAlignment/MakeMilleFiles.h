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
#include <trackbase_historic/ActsTransformations.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrackState_v1.h>

#include <tpc/TpcDistortionCorrectionContainer.h>
#include <tpc/TpcDistortionCorrection.h>
#include <tpc/TpcClusterZCrossingCorrection.h>

class PHCompositeNode;
class PHG4TpcCylinderGeomContainer;
class SvtxTrack;
class SvtxTrackMap;
class TrkrCluster;
class TrkrClusterContainer;
class TpcDistortionCorrectionContainer;
class   ClusterErrorPara;
class Mille;

enum siliconGroup {sensor, stave, barrel};
enum tpcGroup {subsurf, sector, tpc};
enum mmsGroup {tile, mms};

class MakeMilleFiles : public SubsysReco
{
 public:

  MakeMilleFiles(const std::string &name = "MakeMilleFiles");

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void set_datafile_name(const std::string& file) { data_outfilename = file;}
  void set_steeringfile_name(const std::string& file) { steering_outfilename = file;}
  void set_silicon_grouping(int group) {si_group = (siliconGroup) group;}
  void set_tpc_grouping(int group) {tpc_group = (tpcGroup) group;}
  void set_mms_grouping(int group) {mms_group = (mmsGroup) group;}

 private:

Mille* _mille;

int GetNodes(PHCompositeNode* topNode);
Acts::Vector3 getPCALinePoint(Acts::Vector3 global, SvtxTrackState* state);
std::vector<Acts::Vector3> getDerivativesAlignmentAngles(Acts::Vector3& global, TrkrDefs::cluskey cluster_key, TrkrCluster* cluster, Surface surface, int crossing);
SvtxTrack::StateIter getStateIter(Acts::Vector3& global, SvtxTrack* track);
void makeTpcGlobalCorrections(TrkrDefs::cluskey cluster_key, short int crossing, Acts::Vector3& global);
float convertTimeToZ(TrkrDefs::cluskey cluster_key, TrkrCluster *cluster);
Acts::Transform3 makePerturbationTransformation(Acts::Vector3 angles);
int getLabelBase(Acts::GeometryIdentifier id);

  std::map<int, float> derivativeGL;
  std::string  data_outfilename = ("mille_output_data_file.bin");  
  std::string  steering_outfilename = ("steer.txt");  

  /// tpc distortion correction utility class
  TpcDistortionCorrection _distortionCorrection;

  unsigned int _cluster_version = 3;

  ClusterErrorPara _ClusErrPara;

  float sensorAngles[3] = {0.1, 0.1, 0.2};  // perturbation values for each alignment angle

// set default groups to lowest level
  siliconGroup si_group = siliconGroup::sensor;
  tpcGroup tpc_group = tpcGroup::subsurf;
  mmsGroup mms_group = mmsGroup::tile;

  int nstaves[7] = {12,16,20,12,12,16,16};

  std::map<unsigned int, unsigned int> base_layer_map = { {10, 0}, {12,3}, {14,7}, {16,55} };

  SvtxTrackMap *_track_map{nullptr};
  TrkrClusterContainer *_cluster_map{nullptr};						    
  ActsGeometry *_tGeometry{nullptr};

  TpcClusterZCrossingCorrection m_clusterCrossingCorrection;
  TpcDistortionCorrectionContainer* _dcc_static{nullptr};
  TpcDistortionCorrectionContainer* _dcc_average{nullptr};
  TpcDistortionCorrectionContainer* _dcc_fluctuation{nullptr};
};

#endif // MAKEMILLEFILES_H
