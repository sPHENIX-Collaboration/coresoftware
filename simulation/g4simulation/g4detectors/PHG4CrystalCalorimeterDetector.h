#ifndef PHG4CrystalCalorimeterDetector_h
#define PHG4CrystalCalorimeterDetector_h

#include <g4main/PHG4Detector.h>

#include <Geant4/globals.hh>
#include <Geant4/G4Types.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4RotationMatrix.hh>

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

class PHG4CrystalCalorimeterDetector: public PHG4Detector
{

public:

  //! constructor
  PHG4CrystalCalorimeterDetector( PHCompositeNode *Node, const std::string &dnam="BLOCK" );

  //! destructor
  virtual ~PHG4CrystalCalorimeterDetector();

  //! construct
  virtual void Construct( G4LogicalVolume* world );

  //!@name volume accessors
  int IsInCrystalCalorimeter(G4VPhysicalVolume*) const;

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
  int IsActive() const {return _active;}

  int get_Layer() const {return 0;}


private:

  int ConstructCrystals(G4LogicalVolume* envelope);

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
  G4double _crystal_front_dx;
  G4double _crystal_front_dy;
  G4double _crystal_dz;

  G4String _materialCrystal;

  int _active;

  std::string _crystallogicnameprefix;

};

#endif
