// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHMICROMEGASTPCTRACKMATCHING_H
#define PHMICROMEGASTPCTRACKMATCHING_H

#include <fun4all/SubsysReco.h>
#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/ActsSurfaceMaps.h>
#include <tpc/TpcDistortionCorrection.h>
#include <tpc/TpcClusterZCrossingCorrection.h>

#include <array>
#include <string>
#include <vector>

class TrkrClusterContainer;
class TrkrClusterIterationMapv1;
class TrackSeedContainer;
class PHCompositeNode;
class PHG4CylinderGeomContainer;
class TrackSeed;
class TrkrCluster;
class TF1;
class TH1;

class PHMicromegasTpcTrackMatching : public SubsysReco
{
  
  public:
  PHMicromegasTpcTrackMatching(const std::string &name = "PHMicromegasTpcTrackMatching");
  ~PHMicromegasTpcTrackMatching() override = default;

  void set_rphi_search_window_lyr1(const double win){_rphi_search_win[0] = win;}
  void set_z_search_window_lyr1(const double win){_z_search_win[0] = win;}
  void set_rphi_search_window_lyr2(const double win){_rphi_search_win[1] = win;}
  void set_z_search_window_lyr2(const double win){_z_search_win[1] = win;}
  void set_min_tpc_layer(const unsigned int layer){_min_tpc_layer = layer;}
  void set_test_windows_printout(const bool test){_test_windows = test;}
  void set_sc_calib_mode(const bool mode){_sc_calib_mode = mode;}
  void set_collision_rate(const double rate){_collision_rate = rate;}
  void SetIteration(int iter){_n_iteration = iter;}
  void set_track_map_name(const std::string &map_name) { _track_map_name = map_name; }
  
  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;
  
  private:

  //! load nodes relevant for the analysis
  int GetNodes(PHCompositeNode* topNode);

  void copyMicromegasClustersToCorrectedMap( );
  Acts::Vector3 getGlobalPosition( TrkrDefs::cluskey key, TrkrCluster* cluster, short int crossing, unsigned int side);
    
  //! number of layers in the micromegas
  static constexpr unsigned int _n_mm_layers = 2;
  
  bool _use_truth_clusters = false;
  TrkrClusterContainer *_cluster_map{nullptr};
  TrkrClusterContainer *_corrected_cluster_map{nullptr};

  TrackSeedContainer *_track_map{nullptr};

  //! default rphi search window for each layer
  std::array<double,_n_mm_layers> _rphi_search_win = {0.25, 13.0}; 

  //! default z search window for each layer
  std::array<double,_n_mm_layers> _z_search_win = {26.0, 0.25};

  // range of TPC layers to use in projection to micromegas
  unsigned int _min_tpc_layer = 38;
  
  /// first micromegas layer
  /** it is reset in ::Setup using actual micromegas geometry */
  unsigned int _min_mm_layer = 55;

  //! true for initial pass with distorted tracks
  bool _sc_calib_mode = false; 

  //! input rate for phi correction
  double _collision_rate = 50e3;  
  
  //! reference rate for phi correction
  double _reference_collision_rate = 50e3;  

  //! internal event number
  int _event = -1;
  
  //! micomegas geometry
  PHG4CylinderGeomContainer* _geomContainerMicromegas = nullptr;
  TrkrClusterIterationMapv1* _iteration_map = nullptr;
  int _n_iteration = 0;
  std::string _track_map_name = "TpcTrackSeedContainer";

  ActsTrackingGeometry *_tGeometry = nullptr;
  ActsSurfaceMaps *_surfmaps = nullptr;

  /// distortion correction container
  TpcDistortionCorrectionContainer* _dcc = nullptr;
 /// tpc distortion correction utility class
  TpcDistortionCorrection _distortionCorrection;
  TpcClusterZCrossingCorrection _clusterCrossingCorrection;

  //! coarse SC correction function
  TF1 *fdrphi{nullptr};
  double _par0 = -0.36619;
  double _par1 = 0.00375714;

  //! true to printout actual residuals for testing
  bool _test_windows = false;   
  
};

#endif // PHMICROMEGASTPCTRACKMATCHING_H
