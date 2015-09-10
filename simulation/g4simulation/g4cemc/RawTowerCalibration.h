#ifndef RawTowerCalibration_H__
#define RawTowerCalibration_H__

#include <fun4all/SubsysReco.h>
#include <string>

#include <phool/PHTimeServer.h>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeom;

class RawTowerCalibration : public SubsysReco {

 public:
  RawTowerCalibration(const std::string& name="RawTowerCalibration");
  virtual ~RawTowerCalibration(){}

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

#endif /* RawTowerCalibration_H__ */
