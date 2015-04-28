#ifndef __PHG4SVTXDIGITIZER__
#define __PHG4SVTXDIGITIZER__

#include "SvtxHitMap.h"

#include <vector>

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>

class PHG4SvtxDigitizer : public SubsysReco
{
 public:

  PHG4SvtxDigitizer(const char * name = "PHG4SvtxDigitizer");
  ~PHG4SvtxDigitizer(){}
  
  //! module initialization
  int Init(PHCompositeNode *topNode){return 0;}
  
  //! run initialization
  int InitRun(PHCompositeNode *topNode);
  
    //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode) {return 0;}
  
  void set_adc_scale(int layer, unsigned int max_adc, float energy_per_adc) {
    _max_adc.insert(std::make_pair(layer,max_adc));
    _energy_scale.insert(std::make_pair(layer,energy_per_adc));
  }
  
 private:

  void CalculateCylinderCellADCScale(PHCompositeNode *topNode);
  void CalculateLadderCellADCScale(PHCompositeNode *topNode);
  
  void DigitizeCylinderCells(PHCompositeNode *topNode);
  void DigitizeLadderCells(PHCompositeNode *topNode);
  void PrintHits(PHCompositeNode *topNode);
  
  // settings
  std::map<int,unsigned int> _max_adc;
  std::map<int,float> _energy_scale;

  // storage
  SvtxHitMap* _hitmap;
  
  PHTimeServer::timer _timer;   ///< Timer
};

#endif
