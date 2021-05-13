// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHMICROMEGASTPCTRACKMATCHING_H
#define PHMICROMEGASTPCTRACKMATCHING_H

#include <trackreco/PHTrackPropagating.h>

#include <array>
#include <string>
#include <vector>

class PHCompositeNode;
class PHG4CylinderGeomContainer;
class SvtxTrack;
class TrkrCluster;
class TF1;
class TH1;
class TrkrHitSetContainer;

class PHMicromegasTpcTrackMatching : public PHTrackPropagating
{
  
  public:
  PHMicromegasTpcTrackMatching(const std::string &name = "PHMicromegasTpcTrackMatching");
  virtual ~PHMicromegasTpcTrackMatching() = default;

  void set_rphi_search_window_lyr1(const double win){_rphi_search_win[0] = win;}
  void set_z_search_window_lyr1(const double win){_z_search_win[0] = win;}
  void set_rphi_search_window_lyr2(const double win){_rphi_search_win[1] = win;}
  void set_z_search_window_lyr2(const double win){_z_search_win[1] = win;}
  void set_min_tpc_layer(const unsigned int layer){_min_tpc_layer = layer;}
  void set_test_windows_printout(const bool test){_test_windows = test;}
  void set_sc_calib_mode(const bool mode){_sc_calib_mode = mode;}
  void set_collision_rate(const double rate){_collision_rate = rate;}

  protected:
  int Setup(PHCompositeNode* topNode) override;
  int Process() override;
  int End() override;
  
  private:

  int GetNodes(PHCompositeNode* topNode);
    
  // default values, can be replaced from the macro, all in cm
  // rhese correspond to the "baseline" configuration tiles
  static constexpr unsigned int _n_mm_layers = 2;
  std::array<double,_n_mm_layers> _rphi_search_win = {0.25, 13.0}; 
  std::array<double,_n_mm_layers> _z_search_win = {26.0, 0.25};

  // range of TPC layers to use in projection to micromegas
  unsigned int _min_tpc_layer = 38;
  
  /// first micromegas layer
  /** it is reset in ::Setup using actual micromegas geometry */
  unsigned int _min_mm_layer = 55;

  bool _sc_calib_mode = false;  // true for initioal pass with distorted tracks

  double _collision_rate = 50e3;  // input rate for phi correction
  double _reference_collision_rate = 50e3;  // reference rate for phi correction

  int _event = -1;
  
  //! micomegas geometry
  PHG4CylinderGeomContainer* _geomContainerMicromegas = nullptr;

  //! current track
  SvtxTrack *_tracklet_tpc{nullptr};

  TF1 *fdrphi{nullptr};
  double _par0 = -0.36619;
  double _par1 = 0.00375714;

  //!@name evaluation
  //@{
  bool _test_windows = false;   
  std::string _evaluation_rootfile = "PHMicromegasTpcTrackMatching.root";
  std::array<TH1*,_n_mm_layers> _rphi_residuals = {{nullptr, nullptr}};
  std::array<TH1*,_n_mm_layers> _z_residuals = {{nullptr, nullptr}};
  //@}
  
};

#endif // PHMICROMEGASTPCTRACKMATCHING_H
