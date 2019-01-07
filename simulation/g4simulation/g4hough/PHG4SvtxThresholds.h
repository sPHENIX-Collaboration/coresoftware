#ifndef __PHG4SVTXTHRESHOLDS__
#define __PHG4SVTXTHRESHOLDS__

#include <map>
#include <iostream>

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>

class SvtxHitMap;

class PHG4SvtxThresholds : public SubsysReco
{
 public:

  PHG4SvtxThresholds(const std::string &name = "PHG4SvtxThresholds");

  virtual ~PHG4SvtxThresholds(){}
  
  //! module initialization
  int Init(PHCompositeNode *topNode){return 0;}
  
  //! run initialization
  int InitRun(PHCompositeNode *topNode);
  
    //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode);

  void set_threshold(const float fraction_of_mip) {
    std::cout << "PHG4SvtxThresholds use is out of date. "
	      << "Continuing with assumption of tracker with <9 layers. "
	      << "Please update your macros with the latest from GitHub" << std::endl;
    // remove this function eventually
    _fraction_of_mip[0] = fraction_of_mip;
    _fraction_of_mip[1] = fraction_of_mip;
    _fraction_of_mip[2] = fraction_of_mip;
    _fraction_of_mip[3] = fraction_of_mip;
    _fraction_of_mip[4] = fraction_of_mip;
    _fraction_of_mip[5] = fraction_of_mip;
    _fraction_of_mip[6] = fraction_of_mip;
    _fraction_of_mip[7] = fraction_of_mip;
  }
  
  //! set an energy requirement relative to the short-axis MIP expectation
  void set_threshold(const int layer, const float fraction_of_mip) {
    _fraction_of_mip.insert(std::make_pair(layer,fraction_of_mip));
  }
  float get_threshold_by_layer(const int layer) const {
    if (_thresholds_by_layer.find(layer) == _thresholds_by_layer.end()) return 0.0;
    return _thresholds_by_layer.find(layer)->second;
  }

  //! decide if the MIP should use the layer thickness instead of the short-axis
  //! thickness setting in outer layers will kill bending lower pT tracks
  void set_use_thickness_mip(const int layer, const bool use_thickness_mip) {
    _use_thickness_mip.insert(std::make_pair(layer,use_thickness_mip));
  }
  bool get_use_thickness_mip(const int layer) {
    if (_use_thickness_mip.find(layer) == _use_thickness_mip.end()) return false;
    return _use_thickness_mip[layer];
  }

 private:

  void CalculateCylinderThresholds(PHCompositeNode *topNode);
  void CalculateLadderThresholds(PHCompositeNode *topNode);
  void CalculateMVTXLadderThresholds(PHCompositeNode *topNode);

  // settings
  std::map<int,float> _fraction_of_mip;
  std::map<int,float> _thresholds_by_layer;
  std::map<int,bool>  _use_thickness_mip;

  SvtxHitMap* _hits;
  
  PHTimeServer::timer _timer;   ///< Timer
};

#endif
