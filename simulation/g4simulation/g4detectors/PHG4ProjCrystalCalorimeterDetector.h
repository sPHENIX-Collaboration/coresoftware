// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4PROJCRYSTALCALORIMETERDETECTOR_H
#define G4DETECTORS_PHG4PROJCRYSTALCALORIMETERDETECTOR_H

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
  PHG4ProjCrystalCalorimeterDetector( PHG4CrystalCalorimeterSubsystem *subsys, PHCompositeNode *Node, const std::string &dnam="BLOCK" );

  //! destructor
  virtual ~PHG4ProjCrystalCalorimeterDetector();

  //! construct
  virtual void Construct( G4LogicalVolume* world );

  //!@name volume accessors
  virtual int IsInCrystalCalorimeter(G4VPhysicalVolume*) const;

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

  void GetCarbonFiberAdjustments(G4double& adjust_width, G4double& adjust_length);

  void GetCarbonFiberSpacing(G4double& CF_width, G4double& Air_CF, G4double& Air_Cry);
  virtual void SetSupermoduleGeometry(const std::string & filename2)
  {
    _4x4_construct_file = filename2;
  }

private:

  int ConstructProjectiveCrystals(G4LogicalVolume* envelope);
  int Fill4x4Unit(G4LogicalVolume *crystal_logic);
  int FillSpecialUnit(G4LogicalVolume *crystal_logic, G4int ident);


  /* crystal geometry */
  G4double _dx_front;
  G4double _dy_front;
  G4double _dx_back;
  G4double _dy_back;
  G4double _dz_crystal;


  std::string _crystallogicnameprefix;
  std::string _4x4_construct_file;

  bool _overlapcheck_local;
};

#endif
