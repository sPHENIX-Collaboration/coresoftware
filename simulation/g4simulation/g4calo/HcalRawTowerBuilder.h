#ifndef HcalRawTowerBuilder_H__
#define HcalRawTowerBuilder_H__

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <phparameter/PHParameterInterface.h>

#include <string>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeomContainer;

class HcalRawTowerBuilder : public SubsysReco, public PHParameterInterface
{
 public:
  HcalRawTowerBuilder(const std::string &name = "HcalRawTowerBuilder");
  virtual ~HcalRawTowerBuilder() {}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  void Detector(const std::string &d) { detector = d; }
  void EminCut(const double e) { emin = e; }
  void checkenergy(const int i = 1) { chkenergyconservation = i; }

  enum enu_tower_energy_src
  {
    //! save Geant4 energy deposition as the weight of the cells
    kEnergyDeposition,

    //! save light yield as the weight of the cells
    kLightYield,

    //! save ionization energy
    kIonizationEnergy,
    //! initialization value
    unknown = -1
  };

  int get_tower_energy_src() const
  {
    return _tower_energy_src;
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

  void SetDefaultParameters();

 protected:
  void CreateNodes(PHCompositeNode *topNode);
  void ReadParamsFromNodeTree(PHCompositeNode *topNode);

  RawTowerContainer *_towers;
  RawTowerGeomContainer *rawtowergeom;

  std::string detector;
  std::string TowerNodeName;
  std::string TowerGeomNodeName;
  std::string _sim_tower_node_prefix;

  double emin;
  int chkenergyconservation;
  int _tower_energy_src;
  int ncell_to_tower;
  PHTimeServer::timer _timer;
};

#endif /* HcalRawTowerBuilder_H__ */
