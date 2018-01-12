#ifndef Prototype2RawTowerBuilder_H__
#define Prototype2RawTowerBuilder_H__

#include <fun4all/SubsysReco.h>

#include <phool/PHTimeServer.h>

#include <string>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeomContainer;

class Prototype2RawTowerBuilder : public SubsysReco {

 public:
  Prototype2RawTowerBuilder(const std::string& name="Prototype2RawTowerBuilder");
  virtual ~Prototype2RawTowerBuilder(){}

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
    kLightYield,

    //! save ionization energy
    kIonizationEnergy
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

  short get_tower_row(const short cellrow) const;

 protected:
  void CreateNodes(PHCompositeNode *topNode);

  RawTowerContainer* _towers;
  RawTowerGeomContainer *rawtowergeom;

  std::string detector;
  std::string TowerNodeName;
  std::string TowerGeomNodeName;
  std::string _sim_tower_node_prefix;

  double emin;	
  int chkenergyconservation;
  enu_tower_energy_src _tower_energy_src;
  int ncell_to_tower;
  PHTimeServer::timer _timer;

};

#endif /* Prototype2RawTowerBuilder_H__ */
