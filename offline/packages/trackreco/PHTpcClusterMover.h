// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 *  \file		  PHTpcClusterMover.h
 *  \brief		Class for moving corrected TPC clusters to the nearest TPC readout layer radius
 *  \author	 Tony Frawley <afrawley@fsu.edu>
 */

#ifndef PHTPCCLUSTERMOVER_H
#define PHTPCCLUSTERMOVER_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

#include <trackbase/ActsGeometry.h>
#include <trackbase_historic/ActsTransformations.h>
#include <tpc/TpcDistortionCorrectionContainer.h>
#include <tpc/TpcDistortionCorrection.h>

class PHCompositeNode;
class PHG4TpcCylinderGeomContainer;
class SvtxTrack;
class SvtxTrackMap;
class TrkrCluster;
class TrkrClusterContainer;

class PHTpcClusterMover : public SubsysReco
{
 public:

  PHTpcClusterMover(const std::string &name = "PHTpcClusterMover");

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;


 private:

  int GetNodes(PHCompositeNode* topNode);
  int get_circle_circle_intersection(double target_radius, double R, double X0, double Y0, double xref, double yref, double &x, double &y);

 /// acts transformation object
  ActsTransformations _transformer;
  
  /// tpc distortion correction utility class
  TpcDistortionCorrection _distortionCorrection;

  double _z_start=0.0; 
  double _y_start=0.0; 
  double _x_start=0.0; 

  double _z_proj=0.0; 
  double _y_proj=0.0; 
  double _x_proj=0.0; 

  // range of TPC layers to use in projection to micromegas

  PHG4TpcCylinderGeomContainer* _tpc_geom_container = nullptr;

  SvtxTrackMap *_track_map{nullptr};
  SvtxTrack *_track{nullptr};
  TrkrClusterContainer *_cluster_map{nullptr};						    
  TrkrClusterContainer *_corrected_cluster_map{nullptr};						    
  ActsGeometry *_tGeometry{nullptr};
  TpcDistortionCorrectionContainer* _dcc{nullptr};

  double layer_radius[48] = {0};
  double inner_tpc_min_radius = 30.0;
  double mid_tpc_min_radius = 40.0;
  double outer_tpc_min_radius = 60.0;
  double outer_tpc_max_radius = 77.0;

  double inner_tpc_spacing = 0.0;
  double mid_tpc_spacing = 0.0;
  double outer_tpc_spacing = 0.0;

};

#endif // PHTPCCLUSTERMOVER_H
