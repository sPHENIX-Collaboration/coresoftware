// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4CALO_PROTOTYPE2RAWTOWERBUILDER_H
#define G4CALO_PROTOTYPE2RAWTOWERBUILDER_H

#include <phparameter/PHParameterInterface.h>

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class RawTowerGeomContainer;

class Prototype2RawTowerBuilder : public SubsysReco, public PHParameterInterface
{
 public:
  Prototype2RawTowerBuilder(const std::string &name = "Prototype2RawTowerBuilder");
  virtual ~Prototype2RawTowerBuilder() {}
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  void Print(const std::string &what = "ALL") const;

  void SetDefaultParameters();

  void Detector(const std::string &d) { m_Detector = d; }
  void checkenergy(const int i = 1) { m_CheckEnergyConservationFlag = i; }
  enum enu_tower_energy_src
  {
    //! save Geant4 energy deposition as the weight of the cells
    kEnergyDeposition,

    //! save light yield as the weight of the cells
    kLightYield,

    //! save ionization energy
    kIonizationEnergy
  };

  enu_tower_energy_src get_tower_energy_src() const { return m_TowerEnergySrc; }

  void set_tower_energy_src(const enu_tower_energy_src towerEnergySrc)
  {
    m_TowerEnergySrc = towerEnergySrc;
  }

  std::string get_sim_tower_node_prefix() const { return m_SimTowerNodePrefix; }

  void set_sim_tower_node_prefix(const std::string &simTowerNodePrefix)
  {
    m_SimTowerNodePrefix = simTowerNodePrefix;
  }

  short get_tower_row(const short cellrow) const;

 private:
  void ReadParamsFromNodeTree(PHCompositeNode *topNode);

  std::string m_Detector;
  std::string m_TowerNodeName;
  std::string m_TowerGeomNodeName;
  std::string m_SimTowerNodePrefix;

  double m_Emin;
  int m_CheckEnergyConservationFlag;
  enu_tower_energy_src m_TowerEnergySrc;
  int m_NumCellToTower;
};

#endif  // G4CALO_PROTOTYPE2RAWTOWERBUILDER_H
