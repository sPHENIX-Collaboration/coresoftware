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

class PHG4TpcCylinderGeomContainer;

class TpcClusterMover
{
  public:
  
  //! constructor
  TpcClusterMover();

  void set_verbosity(int verb) { _verbosity = verb; }

  std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> processTrack(std::vector<std::pair<TrkrDefs::cluskey,Acts::Vector3>> global_in );

  //! Updates the assumed default geometry below to that contained in the
  //! cell geo
  void initialize_geometry(PHG4TpcCylinderGeomContainer* cellgeo);

  private:
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
  double outer_tpc_max_radius = 76.4;

  double inner_tpc_spacing = 0.0;
  double mid_tpc_spacing = 0.0;
  double outer_tpc_spacing = 0.0;

  int _verbosity = 0;
};

#endif
