#ifndef PHG4ProjCrystalCalorimeterDetector_h
#define PHG4ProjCrystalCalorimeterDetector_h

#include "PHG4CrystalCalorimeterDetector.h"

#include <g4main/PHG4Detector.h>

#include <Geant4/globals.hh>
#include <Geant4/G4Types.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4Material.hh>

#include <string>
#include <map>
#include <vector>
#include <set>

class G4AssemblyVolume;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;

/**
 * \file ${file_name}
 * \brief Module to build crystal calorimeter (endcap) in Geant4
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


  //! Select mapping file for calorimeter tower
  virtual void SetTowerMappingFile( std::string filename ) {
    _inputFile = filename;
  }

  //! Select mapping file for supermodule
  virtual void SetSupermoduleGeometry( std::string filename2 ) {
    _inputFile_4x4_construct = filename2;
  }

  virtual void CrystalDimensions(G4double& dx_front, G4double& dy_front, G4double& dx_back, G4double& dy_back, G4double& dz);

  virtual void SetDimensions(G4double dx_front, G4double dy_front, G4double dx_back, G4double dy_back, G4double dz) {
  _dx_front = dx_front;
  _dy_front = dy_front;
  _dx_back = dx_back;
  _dy_back = dy_back;
  _dz_crystal = dz;
  }

  virtual void CarbonFiberAdjustments(G4double& adjust_width, G4double& adjust_length);

  virtual void CarbonFiberSpacing(G4double& CF_width, G4double& Air_CF, G4double& Air_Cry);

  void SetPlace( G4double place_in_x, G4double place_in_y, G4double place_in_z) {
    _place_in_x = place_in_x;
    _place_in_y = place_in_y;
    _place_in_z = place_in_z;
  }
  void SetXRot( G4double rot_in_x ) { _rot_in_x = rot_in_x; }
  void SetYRot( G4double rot_in_y ) { _rot_in_y = rot_in_y; }
  void SetZRot( G4double rot_in_z ) { _rot_in_z = rot_in_z; }

  void SetMaterial( G4String material ) { _materialCrystal = material; }

  void SetActive(const int i = 1) {_active = i;}
  void SetAbsorberActive(const int i = 1) {_absorberactive = i;}

  void SetInput( G4String inFile ) { _inputFile = inFile; }

  int IsActive() const {return _active;}

  void SuperDetector(const std::string &name) {_superdetector = name;}
  const std::string SuperDetector() const {return _superdetector;}

  int get_Layer() const {return _layer;}

  void BlackHole(const int i=1) {_blackhole = i;}
  int IsBlackHole() const {return _blackhole;}

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
  std::string _inputFile;
  std::string _inputFile_4x4_construct;
};

#endif
