#ifndef RawTowerDigitizer_H__
#define RawTowerDigitizer_H__

#include <fun4all/SubsysReco.h>
#include <string>

#include <phool/PHTimeServer.h>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeom;

class RawTowerDigitizer : public SubsysReco {

 public:
  RawTowerDigitizer(const std::string& name="RawTowerDigitizer");
  virtual ~RawTowerDigitizer(){}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  void Detector(const std::string &d) {detector = d;}
  void EminCut(const double e) {emin = e;}
  void checkenergy(const int i = 1) {chkenergyconservation = i;}

 protected:
  void CreateNodes(PHCompositeNode *topNode);

  RawTowerContainer* _towers;
  RawTowerGeom *rawtowergeom;

  std::string detector;
  std::string TowerNodeName;
  std::string TowerGeomNodeName;

  int _cell_binning;
  double emin;	
  int chkenergyconservation;
  int _nlayers;
  int _nphibins;
  int _netabins;
  double _etamin;
  double _phimin;
  double _etastep;
  double _phistep;

  PHTimeServer::timer _timer;

};

#endif /* RawTowerDigitizer_H__ */
