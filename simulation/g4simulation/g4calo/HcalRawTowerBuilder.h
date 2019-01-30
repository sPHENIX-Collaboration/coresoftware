#ifndef G4CALO_HCALRAWTOWERBUILDER_H
#define G4CALO_HCALRAWTOWERBUILDER_H

#include <fun4all/SubsysReco.h>

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
  void Detector(const std::string &d) { m_Detector = d; }
  void EminCut(const double e) { m_Emin = e; }
  void checkenergy(const int i = 1) { m_ChkEnergyConservationFlag = i; }

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
    return m_TowerEnergySrc;
  }

  std::string
  get_sim_tower_node_prefix() const
  {
    return m_SimTowerNodePrefix;
  }

  void
  set_sim_tower_node_prefix(std::string simTowerNodePrefix)
  {
    m_SimTowerNodePrefix = simTowerNodePrefix;
  }

  short get_tower_row(const short cellrow) const;

  void SetDefaultParameters();

 private:
  void CreateNodes(PHCompositeNode *topNode);
  void ReadParamsFromNodeTree(PHCompositeNode *topNode);

  RawTowerContainer *m_Towers;
  RawTowerGeomContainer *m_RawTowerGeom;

  std::string m_Detector;
  std::string m_TowerNodeName;
  std::string m_TowerGeomNodeName;
  std::string m_SimTowerNodePrefix;

  double m_Emin;
  int m_ChkEnergyConservationFlag;
  int m_TowerEnergySrc;
  int m_NcellToTower;
};

#endif /* G4CALO_HCALRAWTOWERBUILDER_H */
