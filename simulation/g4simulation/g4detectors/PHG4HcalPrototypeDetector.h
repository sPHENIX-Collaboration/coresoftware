// This is the header file for the hcal prototype
// created on 1/27/2014, Liang, HeXC
// Updated on 3/21/2014, Liang, HeXC

#ifndef G4DETECTORS_PHG4HCALPROTOTYPEDETECTOR_H
#define G4DETECTORS_PHG4HCALPROTOTYPEDETECTOR_H

#include "PHG4HcalPrototypeDetectorMessenger.h"
#include <g4main/PHG4Detector.h>

#include <Geant4/globals.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Types.hh>

#include <map>
#include <vector>

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4PVPlacement;
class G4Material;
class G4Box;

class PHG4HcalPrototypeDetector: public PHG4Detector
{

  public:

  //! constructor
  PHG4HcalPrototypeDetector( PHCompositeNode *Node, const std::string &dnam="HCAL", const int lyr = 0 );

  //! destructor
  virtual ~PHG4HcalPrototypeDetector( void )
  {}

  //! construct
  virtual void Construct( G4LogicalVolume* world );

  //!@name volume accessors
  //@{
  int IsInHcalPrototype(G4VPhysicalVolume*) const;
  //@}

  /*
  void SetPlaceZ(const G4double place_z) {place_in_z = place_z*cm;}
  void SetPlace(const G4double place_x, const G4double place_y, const G4double place_z)
  {
    place_in_x = place_x*cm;
    place_in_y = place_y*cm;
    place_in_z = place_z*cm;
  }
  */

  // We will keep these functions for now and deal with them later
  //  void SetXRot(const G4double angle) {x_rot = angle*rad;}
  void SetYRot(const G4double angle) {hcalBoxRotationAngle_z = angle*rad;}
  void SetZRot(const G4double angle) {hcalBoxRotationAngle_y = angle*rad;}
  void SetActive(const int i = 1) {active = i;}
  void SetAbsorberActive(const int i = 1) {absorberactive = i;}
  int IsActive() const {return active;}
  void SuperDetector(const std::string &name) {superdetector = name;}
  const std::string SuperDetector() const {return superdetector;}
  int get_Layer() const {return layer;}

  void BlackHole(const int i=1) {blackhole = i;}
  int IsBlackHole() const {return blackhole;}

  // These functions are copied from our standalone simulation
  void SetMaterial (G4String);            
  void SetOuterHcalDPhi(G4double);
  void SetInnerHcalDPhi(G4double);
  void SetOuterPlateTiltAngle(G4double);
  void SetInnerPlateTiltAngle(G4double);
  void UpdateGeometry();

  private:

  void CalculateGeometry();

  /*
 G4double hcalBoxSizeX, hcalBoxSizeY, hcalBoxSizeZ, hcalBoxRotationAngle_z, hcalBoxRotationAngle_y;

  G4double hcal2SizeX, hcal2SizeY, hcal2SizeZ;
  G4double hcal1SizeX, hcal1SizeY, hcal1SizeZ;

  G4double hcal2RadiusIn, hcal2RadiusOut;
  G4double hcal1RadiusIn, hcal1RadiusOut;

  G4double hcal2Abs_dxa, hcal2Abs_dxb;
  G4double hcal2Abs_dya, hcal2Abs_dyb;
  G4double hcal2Abs_dz;

  G4double hcal1Abs_dxa, hcal1Abs_dxb;
  G4double hcal1Abs_dya, hcal1Abs_dyb;
  G4double hcal1Abs_dz;

  G4int nScint360;
  G4int nHcal2Layers, nHcal1Layers;

  G4double hcal2ScintSizeX, hcal2ScintSizeY, hcal2ScintSizeZ;
  G4double hcal1ScintSizeX, hcal1ScintSizeY, hcal1ScintSizeZ;
  G4double hcal2TiltAngle, hcal1TiltAngle;
  G4double hcal2DPhi, hcal1DPhi;
  G4double hcalJunctionSizeX, hcalJunctionSizeY, hcalJunctionSizeZ;

   */
  G4double hcalBoxSizeX, hcalBoxSizeY, hcalBoxSizeZ, hcalBoxRotationAngle_z, hcalBoxRotationAngle_y;

  G4double hcal2SizeX, hcal2SizeY, hcal2SizeZ;
  G4double hcal1SizeX, hcal1SizeY, hcal1SizeZ;

  G4double hcal2RadiusIn, hcal2RadiusOut;
  G4double hcal1RadiusIn, hcal1RadiusOut;

  G4int nScint360;
  G4int nHcal2Layers, nHcal1Layers;

  G4double hcal2ScintSizeX, hcal2ScintSizeY, hcal2ScintSizeZ;
  G4double hcal1ScintSizeX, hcal1ScintSizeY, hcal1ScintSizeZ;
  G4double hcal2TiltAngle, hcal1TiltAngle;
  G4double hcal2DPhi, hcal1DPhi;

  G4double hcal2Abs_dxa, hcal2Abs_dxb;
  G4double hcal2Abs_dya, hcal2Abs_dyb;
  G4double hcal2Abs_dz;

  G4double hcal1Abs_dxa, hcal1Abs_dxb;
  G4double hcal1Abs_dya, hcal1Abs_dyb;
  G4double hcal1Abs_dz;  

  G4double hcalJunctionSizeX, hcalJunctionSizeY, hcalJunctionSizeZ;

  
  G4VPhysicalVolume* physiWorld;
  G4LogicalVolume*   logicWorld;
  //G4Tubs*            solidWorld;

  G4LogicalVolume*   logicHcalBox;
  G4Box*             solidHcalBox;
  G4PVPlacement*     physiHcalBox;
  

  G4LogicalVolume   *logicHcal2ScintLayer, *logicHcal1ScintLayer;
  G4LogicalVolume   *logicHcal2AbsLayer, *logicHcal1AbsLayer;

  G4Material        *world_mat, *steel, *scint_mat;

  void DefineMaterials();
  G4VPhysicalVolume* ConstructDetector(); 
    
  void SetTiltViaNcross(const int ncross);

  PHG4HcalPrototypeDetectorMessenger*  fDetectorMessenger;

  double outer_tilt_angle;

  int active;
  int absorberactive;
  int layer;
  int blackhole;
  std::string detector_type;
  std::string superdetector;
  
};

#endif
