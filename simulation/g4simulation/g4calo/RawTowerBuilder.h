#ifndef RAWTOWERBUILDER_H__
#define RAWTOWERBUILDER_H__

#include <fun4all/SubsysReco.h>


#include <phool/PHTimeServer.h>
#include <string>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeomContainer;

class RawTowerBuilder : public SubsysReco {

 public:
  RawTowerBuilder(const std::string& name="RawTowerBuilder");
  virtual ~RawTowerBuilder(){}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
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

  enu_tower_energy_src
  get_tower_energy_src() const
  {
    return _tower_energy_src;
  }

  void
  set_tower_energy_src(enu_tower_energy_src towerEnergySrc)
  {
    _tower_energy_src = towerEnergySrc;
  }


  std::string
  get_sim_tower_node_prefix() const
  {
    return _sim_tower_node_prefix;
  }

  void
  set_sim_tower_node_prefix(std::string simTowerNodePrefix)
  {
    _sim_tower_node_prefix = simTowerNodePrefix;
  }

 protected:
  void CreateNodes(PHCompositeNode *topNode);

  RawTowerContainer* _towers;
  RawTowerGeomContainer *rawtowergeom;

  std::string detector;
  std::string TowerNodeName;
  std::string TowerGeomNodeName;
  std::string _sim_tower_node_prefix;

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
