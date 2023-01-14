#ifndef G4CALO_RAWTOWERBUILDER_H
#define G4CALO_RAWTOWERBUILDER_H

#include <fun4all/SubsysReco.h>

#include <g4detectors/PHG4CellDefs.h>

#include <cmath>
#include <string>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeomContainer;

class RawTowerBuilder : public SubsysReco
{
 public:
  RawTowerBuilder(const std::string &name = "RawTowerBuilder");
  ~RawTowerBuilder() override {}
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
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

  RawTowerContainer *m_TowerContainer = nullptr;
  RawTowerGeomContainer *m_RawTowerGeomContainer = nullptr;

  std::string m_Detector = "NONE";
  std::string m_TowerNodeName;
  std::string m_TowerGeomNodeName;
  std::string m_SimTowerNodePrefix;

  enu_tower_energy_src m_TowerEnergySrcEnum = kLightYield;
  int m_CellBinning = PHG4CellDefs::undefined;
  int m_ChkEnergyConservationFlag = 0;
  int m_NumLayers = -1;
  int m_NumPhiBins = -1;
  int m_NumEtaBins = -1;
  double m_Emin = 1e-6;
  double m_EtaMin = NAN;
  double m_PhiMin = NAN;
  double m_EtaStep = NAN;
  double m_PhiStep = NAN;
};

#endif  // G4CALO_RAWTOWERBUILDER_H
