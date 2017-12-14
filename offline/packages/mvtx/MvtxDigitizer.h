#ifndef __MVTXDIGITIZER_H__
#define __MVTXDIGITIZER_H__

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>

#include <vector>

class TrackerHitContainer;

class MvtxDigitizer : public SubsysReco
{
 public:

  MvtxDigitizer(const std::string &name = "MvtxDigitizer");
  virtual ~MvtxDigitizer(){}
  
  //! module initialization
  int Init(PHCompositeNode *topNode){return 0;}
  
  //! run initialization
  int InitRun(PHCompositeNode *topNode);
  
    //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode) {return 0;}
  
  //! set an energy requirement relative to the short-axis MIP expectation
  void set_threshold(const int layer, const float fraction_of_mip) {
    _fraction_of_mip.insert(std::make_pair(layer,fraction_of_mip));
  }
  float get_threshold_by_layer(const int layer) const {
    if (_thresholds_by_layer.find(layer) == _thresholds_by_layer.end()) return 0.0;
    return _thresholds_by_layer.find(layer)->second;
  }

 private:

  void CalculateMapsLadderThresholds(PHCompositeNode *topNode);

  void DigitizeMapsLadderCells(PHCompositeNode *topNode);
  void PrintHits(PHCompositeNode *topNode);
  
  // settings
  std::map<int,float> _fraction_of_mip;
  std::map<int,float> _thresholds_by_layer;

  // storage
  TrackerHitContainer* _hitmap;
  
  PHTimeServer::timer _timer;   ///< Timer
};

#endif
