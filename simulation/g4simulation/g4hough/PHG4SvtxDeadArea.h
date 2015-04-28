#ifndef __PHG4SVTXDEADAREA_H__
#define __PHG4SVTXDEADAREA_H__

#include "SvtxHitMap.h"

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <TRandom3.h>

#include <map>

class PHG4SvtxDeadArea : public SubsysReco {

public:

  PHG4SvtxDeadArea(const char * name = "PHG4SvtxDeadArea");
  ~PHG4SvtxDeadArea(){
    if (_rand) delete _rand;
  }
  
  //! module initialization
  int Init(PHCompositeNode *topNode){return 0;}
  
  //! run initialization
  int InitRun(PHCompositeNode *topNode);
  
    //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode);
  
  //! kill random hits in this layer with fractional eff efficiency
  void set_hit_efficiency(int ilayer, float eff) {
    _eff_by_layer.insert(std::make_pair(ilayer,eff));
  }
  float get_hit_efficiency(int ilayer) {
    if (_eff_by_layer.find(ilayer) == _eff_by_layer.end()) return 1.0;
    return _eff_by_layer[ilayer];
  }

 private:

  void FillCylinderDeadAreaMap(PHCompositeNode *topNode);
  void FillLadderDeadAreaMap(PHCompositeNode *topNode);

  // settings
  std::map<int,float> _eff_by_layer;

  // storage
  SvtxHitMap* _hits;
  TRandom3* _rand;

  PHTimeServer::timer _timer;   ///< Timer
};

#endif
