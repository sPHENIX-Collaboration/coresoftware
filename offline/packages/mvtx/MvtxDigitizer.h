#ifndef __MVTXDIGITIZER_H__
#define __MVTXDIGITIZER_H__

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <tracker/TrkrHitSetContainer.h>

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
  void SetThreshold(const int layer, const float fraction_of_mip) {
    fraction_of_mip_.insert(std::make_pair(layer,fraction_of_mip));
  }
  float GetThresholdByLayer(const int layer) const {
    if (thresholds_by_layer_.find(layer) == thresholds_by_layer_.end()) return 0.0;
    return thresholds_by_layer_.find(layer)->second;
  }

 private:

  void CalculateMapsLadderThresholds(PHCompositeNode *topNode);

  void DigitizeMapsLadderCells(PHCompositeNode *topNode);
  void PrintHits(PHCompositeNode *topNode);
  
  // settings
  std::map<int,float> fraction_of_mip_;
  std::map<int,float> thresholds_by_layer_;

  // storage
  TrkrHitSetContainer* hitmap_;
  
  PHTimeServer::timer timer_;   ///< Timer
};

#endif
