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

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }

  int get_Layer() const { return _layer; }

  // ----- additional accessors used by derived classes: ------------

  //! Select mapping file for supermodule
  virtual void SetSupermoduleGeometry(const std::string &filename2)
  {
    return;
  }

  PHParameters *GetParams() {return m_Params;}

 protected:  // for variable also used in PHG4ProjCrystalCalorimeterDetector
  PHG4CrystalCalorimeterDisplayAction *GetDisplayAction() { return m_DisplayAction; }
  /* Calorimeter envelope geometry */

  G4double _rMin1;
  G4double _rMax1;
  G4double _rMin2;
  G4double _rMax2;

  G4double _dZ;
  G4double _sPhi;
  G4double _dPhi;

  /* general detector parameters */
  int _layer;

  std::string m_SuperDetector;
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

  PHG4CrystalCalorimeterDisplayAction *m_DisplayAction = nullptr;

  std::string _towerlogicnameprefix;

  std::map<std::string, G4double> _map_global_parameter;
  std::map<std::string, towerposition> _map_tower;
  std::set<G4VPhysicalVolume *> m_ActiveVolumeSet;
  std::set<G4VPhysicalVolume *> m_PassiveVolumeSet;
  // since getting parameters is a map search we do not want to
  // do this in every step, the parameters used are cached
  // in the following variables
  int m_IsActive;
  int m_AbsorberActive;
};

#endif
