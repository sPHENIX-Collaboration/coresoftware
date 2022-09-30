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


 private:

  Mille *_mille;

  int GetNodes(PHCompositeNode* topNode);
  Acts::Vector3 getPCALinePoint(Acts::Vector3 global, SvtxTrackState* state);

 /// acts transformation object
  //  ActsTransformations _transformer;

  std::map<int, float> derivativeGL;
  
  /// tpc distortion correction utility class
  TpcDistortionCorrection _distortionCorrection;

  //  PHG4TpcCylinderGeomContainer* _tpc_geom_container = nullptr;

  unsigned int _cluster_version = 3;

  ClusterErrorPara _ClusErrPara;

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
