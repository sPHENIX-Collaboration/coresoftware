// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4FORWARDHCALDETECTOR_H
#define G4DETECTORS_PHG4FORWARDHCALDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <map>
#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHG4ForwardHcalDisplayAction;
class PHG4ForwardHcalSubsystem;

/**
 * \file ${file_name}
 * \brief Module to build forward sampling Hadron calorimeterr (endcap) in Geant4
 * \author Nils Feege <nils.feege@stonybrook.edu>
 */

class PHG4ForwardHcalDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4ForwardHcalDetector(PHG4ForwardHcalSubsystem *subsys, PHCompositeNode *Node, const std::string &dnam = "BLOCK");

  //! destructor
  virtual ~PHG4ForwardHcalDetector() {}

  //! construct
  virtual void Construct(G4LogicalVolume *world);

  //!@name volume accessors
  int IsInForwardHcal(G4VPhysicalVolume *) const;

  //! Select mapping file for calorimeter tower
  void SetTowerMappingFile(const std::string &filename)
  {
    _mapping_tower_file = filename;
  }

  void SetTowerDimensions(G4double dx, G4double dy, G4double dz)
  {
    _tower_dx = dx;
    _tower_dy = dy;
    _tower_dz = dz;
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

  void SetMaterialScintillator(G4String material) { _materialScintillator = material; }
  void SetMaterialAbsorber(G4String material) { _materialAbsorber = material; }

  void SetActive(const int i = 1) { _active = i; }
  void SetAbsorberActive(const int i = 1) { _absorberactive = i; }

  int IsActive() const { return _active; }

  void SuperDetector(const std::string &name) { _superdetector = name; }
  const std::string SuperDetector() const { return _superdetector; }

  int get_Layer() const { return _layer; }

  void BlackHole(const int i = 1) { _blackhole = i; }
  int IsBlackHole() const { return _blackhole; }

 private:
  G4LogicalVolume *ConstructTower();
  int PlaceTower(G4LogicalVolume *envelope, G4LogicalVolume *tower);
  int ParseParametersFromTable();

  struct towerposition
  {
    G4double x;
    G4double y;
    G4double z;
  };

  PHG4ForwardHcalDisplayAction *m_DisplayAction;

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

  /* HCAL tower geometry */
  G4double _tower_dx;
  G4double _tower_dy;
  G4double _tower_dz;

  G4String _materialScintillator;
  G4String _materialAbsorber;

  int _active;
  int _absorberactive;
  int _layer;
  int _blackhole;

  std::string _towerlogicnameprefix;
  std::string _superdetector;
  std::string _mapping_tower_file;

  std::map<std::string, G4double> _map_global_parameter;
  std::map<std::string, towerposition> _map_tower;
};

#endif
