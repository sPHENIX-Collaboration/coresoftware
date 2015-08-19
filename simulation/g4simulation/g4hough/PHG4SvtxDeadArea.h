#ifndef __PHG4SVTXDEADAREA_H__
#define __PHG4SVTXDEADAREA_H__

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>

#include <map>

class SvtxHitMap;
class TRandom3;

class PHG4SvtxDeadArea : public SubsysReco {

public:

  PHG4SvtxDeadArea(const std::string &name = "PHG4SvtxDeadArea");
  virtual ~PHG4SvtxDeadArea();
  
  //! module initialization
  int Init(PHCompositeNode *topNode){return 0;}
  
  //! run initialization
  int InitRun(PHCompositeNode *topNode);
  
    //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode);
  
  //! kill random hits in this layer with fractional eff efficiency
  void set_hit_efficiency(const int ilayer, const float eff) {
    _eff_by_layer.insert(std::make_pair(ilayer,eff));
  }
  float get_hit_efficiency(const int ilayer) const {
    if (_eff_by_layer.find(ilayer) == _eff_by_layer.end()) return 1.0;
    return _eff_by_layer.find(ilayer)->second;
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
