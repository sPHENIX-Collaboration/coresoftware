#ifndef RAWTOWERBUILDER_H__
#define RAWTOWERBUILDER_H__

#include <fun4all/SubsysReco.h>
#include <string>

#include <phool/PHTimeServer.h>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeom;

class RawTowerBuilder : public SubsysReco {

 public:
  RawTowerBuilder(const std::string& name="RawTowerBuilder");
  virtual ~RawTowerBuilder(){}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  void Detector(const std::string &d) {detector = d;}
  void EminCut(const double e) {emin = e;}
  void checkenergy(const int i = 1) {chkenergyconservation = i;}

  enum enu_tower_energy_src
  {
    //! save Geant4 energy deposition as the weight of the cells
      kEnergyDeposition,

      //! save light yield as the weight of the cells
      kLightYield

  };


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
  enu_tower_energy_src _tower_energy_src;

  PHTimeServer::timer _timer;

};

#endif /* RAWTOWERBUILDER_H__ */
