// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CRYSTALCALORIMETERDETECTOR_H
#define G4DETECTORS_PHG4CRYSTALCALORIMETERDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4Types.hh>

#include <map>
#include <set>
#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4CrystalCalorimeterDisplayAction;
class PHG4Subsystem;
class PHParameters;

/**
 * \file ${file_name}
 * \brief Module to build crystal calorimeter in Geant4
 * \author Nils Feege <nils.feege@stonybrook.edu>
 */

class PHG4CrystalCalorimeterDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4CrystalCalorimeterDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  virtual ~PHG4CrystalCalorimeterDetector() {}

  //! construct
  virtual void ConstructMe(G4LogicalVolume *world);

  //! check if volume is in this calorimeter
  virtual int IsInCrystalCalorimeter(G4VPhysicalVolume *) const;

  // ----- accessing member variables: ------------

  //! Select mapping file for calorimeter tower
  void SetTowerMappingFile(std::string &filename)
  {
    _mapping_tower_file = filename;
  }

  void SetPlace(G4double place_in_x, G4double place_in_y, G4double place_in_z)
  {
    _place_in_x = place_in_x;
    _place_in_y = place_in_y;
    _place_in_z = place_in_z;
  }

  void SetRotation(G4double rot_in_x, G4double rot_in_y, G4double rot_in_z)
  {
    _rot_in_x = rot_in_x;
    _rot_in_y = rot_in_y;
    _rot_in_z = rot_in_z;
  }

  void SetCrystalSize(G4double dx, G4double dy, G4double dz)
  {
    _crystal_dx = dx;
    _crystal_dy = dy;
    _crystal_dz = dz;
  }

  void SetMaterialCrystal(G4String material)
  {
    _materialCrystal = material;
  }

  void SetActive(const int i = 1) { _active = i; }
  void SetAbsorberActive(const int i = 1) { _absorberactive = i; }

  int IsActive() const { return _active; }

  void SuperDetector(const std::string &name) { _superdetector = name; }
  const std::string SuperDetector() const { return _superdetector; }

  int get_Layer() const { return _layer; }

  void BlackHole(const int i = 1) { _blackhole = i; }
  int IsBlackHole() const { return _blackhole; }

  // ----- additional accessors used by derived classes: ------------

  //! Select mapping file for supermodule
  virtual void SetSupermoduleGeometry(const std::string &filename2)
  {
    return;
  }

 protected:  // for variable also used in PHG4ProjCrystalCalorimeterDetector
  PHG4CrystalCalorimeterDisplayAction *GetDisplayAction() { return m_DisplayAction; }
  /* Calorimeter envelope geometry */
  G4double _place_in_x;
  G4double _place_in_y;
  G4double _place_in_z;

  G4double _rot_in_x;
  G4double _rot_in_y;
  G4double _rot_in_z;

  G4double _rMin1;
  G4double _rMax1;
  G4double _rMin2;
  G4double _rMax2;

  G4double _dZ;
  G4double _sPhi;
  G4double _dPhi;

  /* crystal geometry */
  G4double _crystal_dx;
  G4double _crystal_dy;
  G4double _crystal_dz;

  G4String _materialCrystal;

  /* general detector parameters */
  int _active;
  int _absorberactive;
  int _layer;
  int _blackhole;

  std::string _superdetector;
  std::string _mapping_tower_file;

 private:  // private stuff
  G4LogicalVolume *ConstructTower();
  int PlaceTower(G4LogicalVolume *envelope, G4LogicalVolume *tower);
  int ParseParametersFromTable();

  struct towerposition
  {
    G4double x;
    G4double y;
    G4double z;
  };
  PHParameters *m_Params = nullptr;

  PHG4CrystalCalorimeterDisplayAction *m_DisplayAction;

  std::string _towerlogicnameprefix;

  std::map<std::string, G4double> _map_global_parameter;
  std::map<std::string, towerposition> _map_tower;
  std::set<G4VPhysicalVolume*> m_ActiveVolumeSet;
  std::set<G4VPhysicalVolume*> m_PassiveVolumeSet;
};

#endif
