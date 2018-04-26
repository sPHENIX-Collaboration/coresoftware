// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4CALO_PROTOTYPE2RAWTOWERBUILDER_H
#define G4CALO_PROTOTYPE2RAWTOWERBUILDER_H

#include <g4detectors/PHG4ParameterInterface.h>

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class RawTowerGeomContainer;

class Prototype2RawTowerBuilder : public SubsysReco, public PHG4ParameterInterface
{

 public:
  Prototype2RawTowerBuilder(const std::string& name="Prototype2RawTowerBuilder");
  virtual ~Prototype2RawTowerBuilder(){}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void SetDefaultParameters();

  void Detector(const std::string &d) {m_Detector = d;}
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

  RawTowerGeomContainer *rawtowergeom;

  std::string m_Detector;
  std::string TowerNodeName;
  std::string TowerGeomNodeName;
  std::string _sim_tower_node_prefix;

  double emin;	
  int chkenergyconservation;
  enu_tower_energy_src _tower_energy_src;
  int ncell_to_tower;

};

#endif // G4CALO_PROTOTYPE2RAWTOWERBUILDER_H
