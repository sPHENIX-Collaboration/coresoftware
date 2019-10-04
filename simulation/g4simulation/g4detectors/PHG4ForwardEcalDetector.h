// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4FORWARDECALDETECTOR_H
#define G4DETECTORS_PHG4FORWARDECALDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Types.hh>  // for G4double

#include <map>
#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4ForwardEcalDisplayAction;
class PHG4Subsystem;
class PHG4GDMLConfig;

/**
 * \file ${file_name}
 * \brief Module to build forward sampling Hadron calorimeterr (endcap) in Geant4
 * \author Nils Feege <nils.feege@stonybrook.edu>
 */

class PHG4ForwardEcalDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4ForwardEcalDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, const std::string &dnam);

  //! destructor
  virtual ~PHG4ForwardEcalDetector(){}

  //! construct
  virtual void ConstructMe(G4LogicalVolume *world);

  //!@name volume accessors
  int IsInForwardEcal(G4VPhysicalVolume *) const;

  //! Select mapping file for calorimeter tower
  void SetTowerMappingFile(const std::string &filename)
  {
    _mapping_tower_file = filename;
  }

  void SetTowerDimensions(G4double dx, G4double dy, G4double dz, G4double type)
  {
    if (type == 0)
    {
      _tower0_dx = dx;
      _tower0_dy = dy;
      _tower0_dz = dz;
    }
    else if (type == 1)
    {
      _tower1_dx = dx;
      _tower1_dy = dy;
      _tower1_dz = dz;
    }
    else if (type == 2)
    {
      _tower2_dx = dx;
      _tower2_dy = dy;
      _tower2_dz = dz;
    }
    else if (type == 3)
    {
      _tower3_dx = dx;
      _tower3_dy = dy;
      _tower3_dz = dz;
    }
    else if (type == 4)
    {
      _tower4_dx = dx;
      _tower4_dy = dy;
      _tower4_dz = dz;
    }
    else if (type == 5)
    {
      _tower5_dx = dx;
      _tower5_dy = dy;
      _tower5_dz = dz;
    }
    else if (type == 6)
    {
      _tower6_dx = dx;
      _tower6_dy = dy;
      _tower6_dz = dz;
    }
  }

  void SetPlace(G4double place_in_x, G4double place_in_y, G4double place_in_z)
  {
    _place_in_x = place_in_x;
    _place_in_y = place_in_y;
    _place_in_z = place_in_z;
  }

  void SetXRot(G4double rot_in_x) { _rot_in_x = rot_in_x; }
  void SetYRot(G4double rot_in_y) { _rot_in_y = rot_in_y; }
  void SetZRot(G4double rot_in_z) { _rot_in_z = rot_in_z; }

  void SetActive(const int i = 1) { _active = i; }
  void SetAbsorberActive(const int i = 1) { _absorberactive = i; }

  int IsActive() const { return _active; }

  void SuperDetector(const std::string &name) { _superdetector = name; }
  const std::string SuperDetector() const { return _superdetector; }

  int get_Layer() const { return _layer; }

  void BlackHole(const int i = 1) { _blackhole = i; }
  int IsBlackHole() const { return _blackhole; }
  PHG4ForwardEcalDisplayAction *GetDisplayAction() { return m_DisplayAction; }

 private:
  G4LogicalVolume *ConstructTower(int type);
  G4LogicalVolume *ConstructTowerType2();
  G4LogicalVolume *ConstructTowerType3_4_5_6(int type);
  int PlaceTower(G4LogicalVolume *envelope, G4LogicalVolume *tower[6]);
  int ParseParametersFromTable();

  struct towerposition
  {
    G4double x;
    G4double y;
    G4double z;
    G4double type;
  };

  PHG4ForwardEcalDisplayAction *m_DisplayAction;

  std::map<std::string, towerposition> _map_tower;

  /* ECAL tower geometry */
  G4double _tower0_dx;
  G4double _tower0_dy;
  G4double _tower0_dz;

  G4double _tower1_dx;
  G4double _tower1_dy;
  G4double _tower1_dz;

  G4double _tower2_dx;
  G4double _tower2_dy;
  G4double _tower2_dz;

  G4double _tower3_dx;
  G4double _tower3_dy;
  G4double _tower3_dz;

  G4double _tower4_dx;
  G4double _tower4_dy;
  G4double _tower4_dz;

  G4double _tower5_dx;
  G4double _tower5_dy;
  G4double _tower5_dz;

  G4double _tower6_dx;
  G4double _tower6_dy;
  G4double _tower6_dz;

 protected:
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

  int _active;
  int _absorberactive;
  int _layer;
  int _blackhole;

  std::string _towerlogicnameprefix;
  std::string _superdetector;
  std::string _mapping_tower_file;

  std::map<std::string, G4double> _map_global_parameter;

  //! registry for volumes that should not be exported, i.e. fibers
  PHG4GDMLConfig *gdml_config;
};

#endif
