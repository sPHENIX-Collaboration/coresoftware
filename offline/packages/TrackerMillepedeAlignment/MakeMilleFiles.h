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
class Mille;
class   ClusterErrorPara;

class MakeMilleFiles : public SubsysReco
{
 public:

  MakeMilleFiles(const std::string &name = "MakeMilleFiles");

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void set_datafile_name(std::string file) { data_outfilename = file;}
  void set_steeringfile_name(std::string file) { steering_outfilename = file;}

 private:

  Mille *_mille;

  int GetNodes(PHCompositeNode* topNode);
  Acts::Vector3 getPCALinePoint(Acts::Vector3 global, SvtxTrackState* state);
std::vector<Acts::Vector3> getDerivativesAlignmentAngles(Acts::Vector3& global, TrkrDefs::cluskey cluster_key, TrkrCluster* cluster, Surface surface, int crossing);
  SvtxTrack::StateIter getStateIter(Acts::Vector3& global, SvtxTrack* track);
  void makeTpcGlobalCorrections(TrkrDefs::cluskey cluster_key, short int crossing, Acts::Vector3& global);
  float convertTimeToZ(TrkrDefs::cluskey cluster_key, TrkrCluster *cluster);
  Acts::Transform3 makePerturbationTransformation(Acts::Vector3 angles);

  std::map<int, float> derivativeGL;
  std::string  data_outfilename = ("mille_output_data_file.bin");  
  std::string  steering_outfilename = ("steer.txt");  

  /// tpc distortion correction utility class
  TpcDistortionCorrection _distortionCorrection;

  unsigned int _cluster_version = 3;

  ClusterErrorPara _ClusErrPara;

float sensorAngles[3] = {0.1, 0.1, 0.2};  // perturbation values for each alignment angle


  SvtxTrackMap *_track_map{nullptr};
  SvtxTrack *_track{nullptr};
  TrkrClusterContainer *_cluster_map{nullptr};						    
  ActsGeometry *_tGeometry{nullptr};

  TpcClusterZCrossingCorrection m_clusterCrossingCorrection;
  TpcDistortionCorrectionContainer* _dcc_static{nullptr};
  TpcDistortionCorrectionContainer* _dcc_average{nullptr};
  TpcDistortionCorrectionContainer* _dcc_fluctuation{nullptr};
};

#endif // MAKEMILLEFILES_H
