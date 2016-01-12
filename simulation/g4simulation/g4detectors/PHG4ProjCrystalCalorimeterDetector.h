#ifndef PHG4ProjCrystalCalorimeterDetector_h
#define PHG4ProjCrystalCalorimeterDetector_h

#include "PHG4CrystalCalorimeterDetector.h"

/**
 * \file ${file_name}
 * \brief Module to build crystal calorimeter with projective crystal geometry in Geant4
 * \author Nils Feege <nils.feege@stonybrook.edu>
 */

class PHG4ProjCrystalCalorimeterDetector: public PHG4CrystalCalorimeterDetector
{

public:

  //! constructor
  PHG4ProjCrystalCalorimeterDetector( PHCompositeNode *Node, const std::string &dnam="BLOCK" );

  //! destructor
  virtual ~PHG4ProjCrystalCalorimeterDetector();

  //! construct
  virtual void Construct( G4LogicalVolume* world );

  //!@name volume accessors
  int IsInCrystalCalorimeter(G4VPhysicalVolume*) const;


  // ----- accessing member variables: ------------

  void SetCrystalSize(G4double dx_front, G4double dy_front, G4double dx_back, G4double dy_back, G4double dz) {
  _dx_front = dx_front;
  _dy_front = dy_front;
  _dx_back = dx_back;
  _dy_back = dy_back;
  _dz_crystal = dz;
  }

  void GetCrystalSize(G4double& dx_front, G4double& dy_front, G4double& dx_back, G4double& dy_back, G4double& dz) {
  dx_front = _dx_front;
  dy_front = _dy_front;
  dx_back = _dx_back;
  dy_back = _dy_back;
  dz = _dz_crystal;
  }

  void CarbonFiberAdjustments(G4double& adjust_width, G4double& adjust_length);

  void CarbonFiberSpacing(G4double& CF_width, G4double& Air_CF, G4double& Air_Cry);


private:

  int ConstructProjectiveCrystals(G4LogicalVolume* envelope);
  int Fill4x4Unit(G4LogicalVolume *crystal_logic);
  int FillSpecialUnit(G4LogicalVolume *crystal_logic, G4int ident);

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

  G4String _materialCrystal;

  /* crystal geometry */
  G4double _dx_front;
  G4double _dy_front;
  G4double _dx_back;
  G4double _dy_back;
  G4double _dz_crystal;

  int _active;
  int _absorberactive;
  int _layer;
  int _blackhole;

  std::string _superdetector;

  std::string _crystallogicnameprefix;
  std::string _mapping_tower_file;
  std::string _4x4_construct_file;
};

#endif
