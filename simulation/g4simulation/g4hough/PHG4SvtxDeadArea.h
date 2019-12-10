#ifndef G4HOUGH_PHG4SVTXDEADAREA_H
#define G4HOUGH_PHG4SVTXDEADAREA_H

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>

// rootcint barfs with this header so we need to hide it
#if !defined (__CINT__) || defined(__CLING__)
#include <gsl/gsl_rng.h>
#endif

#include <map>
#include <string>                // for string
#include <utility>               // for pair, make_pair

class PHCompositeNode;
class SvtxHitMap;

class PHG4SvtxDeadArea : public SubsysReco {

public:

  PHG4SvtxDeadArea(const std::string &name = "PHG4SvtxDeadArea");
  virtual ~PHG4SvtxDeadArea();
  
  //! run initialization
  int InitRun(PHCompositeNode *topNode);
  
    //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! kill random hits in this layer with fractional eff efficiency
  void set_hit_efficiency(const int ilayer, const float eff) {
    _eff_by_layer.insert(std::make_pair(ilayer,eff));
  }
  float get_hit_efficiency(const int ilayer) const {
    if (_eff_by_layer.find(ilayer) == _eff_by_layer.end()) return 1.0;
    return _eff_by_layer.find(ilayer)->second;
  }

 private:

  void GenericFillDeadAreaMap(PHCompositeNode *topNode, const std::string &detectorname);

  // settings
  std::map<int,float> _eff_by_layer;

  // storage
  SvtxHitMap* _hits;

  PHTimeServer::timer _timer;   ///< Timer

#ifndef __CINT__
  gsl_rng *RandomGenerator;
#endif

};

#endif
