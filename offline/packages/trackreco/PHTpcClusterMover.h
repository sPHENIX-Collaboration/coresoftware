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

#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase_historic/ActsTransformations.h>
#include <tpc/TpcDistortionCorrectionContainer.h>
#include <tpc/TpcDistortionCorrection.h>

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class TrkrCluster;
class TrkrClusterContainer;

class PHTpcClusterMover : public SubsysReco
{
 public:

  PHTpcClusterMover(const std::string &name = "PHTpcClusterMover");

  ~PHTpcClusterMover() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;


 private:

  int GetNodes(PHCompositeNode* topNode);

void CircleFitByTaubin (std::vector<Acts::Vector3> clusters, double &R, double &X0, double &Y0);
  void circle_circle_intersection(double r1, double r2, double x2, double y2, double &xplus, double &yplus, double &xminus, double &yminus);
void  line_fit(std::vector<Acts::Vector3> clusters, double &a, double &b);
  int get_circle_circle_intersection(double target_radius, double R, double X0, double Y0, double xref, double yref, double &x, double &y);

  Surface get_tpc_surface_from_coords(TrkrDefs::hitsetkey hitsetkey,
					Acts::Vector3 world,
					ActsSurfaceMaps *surfMaps,
					ActsTrackingGeometry *tGeometry,
					TrkrDefs::subsurfkey& subsurfkey);


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

  SvtxTrackMap *_track_map{nullptr};
  SvtxTrack *_track{nullptr};
  TrkrClusterContainer *_cluster_map{nullptr};						    
  TrkrClusterContainer *_corrected_cluster_map{nullptr};						    
  ActsSurfaceMaps *_surfmaps{nullptr};
  ActsTrackingGeometry *_tGeometry{nullptr};
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
