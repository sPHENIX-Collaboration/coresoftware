#ifndef __PHG4SVTXTHRESHOLDS__
#define __PHG4SVTXTHRESHOLDS__

#include <vector>

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
  
  //! set an energy requirement relative to the short-axis MIP expectation
  void set_threshold(const float fraction_of_mip) {
    _fraction_of_mip = fraction_of_mip;
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

  // settings
  float _fraction_of_mip;
  std::map<int,float> _thresholds_by_layer;
  std::map<int,bool> _use_thickness_mip;

  SvtxHitMap* _hits;
  
  PHTimeServer::timer _timer;   ///< Timer
};

#endif
