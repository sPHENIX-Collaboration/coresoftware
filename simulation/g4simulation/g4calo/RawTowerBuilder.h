#ifndef G4CALO_RAWTOWERBUILDER_H
#define G4CALO_RAWTOWERBUILDER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeomContainer;

class RawTowerBuilder : public SubsysReco
{
 public:
  RawTowerBuilder(const std::string &name = "RawTowerBuilder");
  virtual ~RawTowerBuilder() {}
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  void Detector(const std::string &d) { m_Detector = d; }
  void EminCut(const double e) { m_Emin = e; }
  void checkenergy(const int i = 1) { m_ChkEnergyConservationFlag = i; }
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
    return m_TowerEnergySrcEnum;
  }

  void
  set_tower_energy_src(enu_tower_energy_src towerEnergySrc)
  {
    m_TowerEnergySrcEnum = towerEnergySrc;
  }

  std::string
  get_sim_tower_node_prefix() const
  {
    return m_SimTowerNodePrefix;
  }

  void
  set_sim_tower_node_prefix(const std::string &simTowerNodePrefix)
  {
    m_SimTowerNodePrefix = simTowerNodePrefix;
  }

 protected:
  void CreateNodes(PHCompositeNode *topNode);

  RawTowerContainer *m_TowerContainer;
  RawTowerGeomContainer *m_RawTowerGeomContainer;

  std::string m_Detector;
  std::string m_TowerNodeName;
  std::string m_TowerGeomNodeName;
  std::string m_SimTowerNodePrefix;

  enu_tower_energy_src m_TowerEnergySrcEnum;
  int m_CellBinning;
  int m_ChkEnergyConservationFlag;
  int m_NumLayers;
  int m_NumPhiBins;
  int m_NumEtaBins;
  double m_Emin;
  double m_EtaMin;
  double m_PhiMin;
  double m_EtaStep;
  double m_PhiStep;
};

#endif  // G4CALO_RAWTOWERBUILDER_H
