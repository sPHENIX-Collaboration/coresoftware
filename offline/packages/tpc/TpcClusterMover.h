#ifndef TPC_TPCCLUSTERMOVER_H
#define TPC_TPCCLUSTERMOVER_H

/*!
 * \file TpcClusterMover.h
 * \Moves TPC clusters to the readout TPC surface after distortion corrections
 * \author Tony Frawley, May 2022
 */
#include <vector>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/ActsTransformations.h>

class TpcClusterMover
{
  public:
  
  //! constructor
  TpcClusterMover();

  void set_verbosity(int verb) { _verbosity = verb; }
  std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> processTrack(std::vector<std::pair<TrkrDefs::cluskey,Acts::Vector3>> global_in );
  void CircleFitByTaubin (std::vector<Acts::Vector3> clusters, double &R, double &X0, double &Y0);
  void  line_fit(std::vector<Acts::Vector3> clusters, double &a, double &b);
 void circle_circle_intersection(double r1, double r2, double x2, double y2, double &xplus, double &yplus, double &xminus, double &yminus);
 int get_circle_circle_intersection(double target_radius, double R, double X0, double Y0, double xclus, double yclus, double &x, double &y);

 double _z_start=0.0; 
  double _y_start=0.0; 
  double _x_start=0.0; 

  double _z_proj=0.0; 
  double _y_proj=0.0; 
  double _x_proj=0.0; 

  double layer_radius[48] = {0};
  double inner_tpc_min_radius = 30.0;
  double mid_tpc_min_radius = 40.0;
  double outer_tpc_min_radius = 60.0;
  double outer_tpc_max_radius = 77.0;

  double inner_tpc_spacing = 0.0;
  double mid_tpc_spacing = 0.0;
  double outer_tpc_spacing = 0.0;

  int _verbosity = 0;
};

#endif
