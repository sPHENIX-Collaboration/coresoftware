// This is the header file for the 2nd hcal prototype
// created on 11/20/2015, HeXC

#ifndef PHG4HcalPrototype2Detector_h
#define PHG4HcalPrototype2Detector_h

#include "g4main/PHG4Detector.h"
#include "PHG4HcalPrototype2DetectorMessenger.h"

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

class PHG4HcalPrototype2Detector: public PHG4Detector
{

  public:

  //! constructor
  PHG4HcalPrototype2Detector( PHCompositeNode *Node, const std::string &dnam="HCAL", const int lyr = 0 );

  //! destructor
  virtual ~PHG4HcalPrototype2Detector( void )
  {}

  //! construct
  virtual void Construct( G4LogicalVolume* world );

  //!@name volume accessors
  //@{
  int IsInHcalPrototype2(G4VPhysicalVolume*) const;
  //@}

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

  G4double hcal2AbsEx_dxa, hcal2AbsEx_dxb;
  G4double hcal2AbsEx_dya, hcal2AbsEx_dyb;
  G4double hcal2AbsEx_dz;

  G4double hcal1Abs_dxa, hcal1Abs_dxb;
  G4double hcal1Abs_dya, hcal1Abs_dyb;
  G4double hcal1Abs_dz;  

  G4double cryostatSizeX, cryostatSizeY, cryostatSizeZ;
  
  G4VPhysicalVolume* physiWorld;
  G4LogicalVolume*   logicWorld;

  G4LogicalVolume*   logicHcalBox;
  G4Box*             solidHcalBox;
  G4PVPlacement*     physiHcalBox;
  
  G4LogicalVolume   *logicHcal2ScintLayer, *logicHcal1ScintLayer;
  G4LogicalVolume   *logicHcal2AbsLayer, *logicHcal1AbsLayer;
  G4LogicalVolume   *logicHcal1AbsExtra, *logicHcal2AbsExtra;
  G4LogicalVolume   *logicHcal2AbsExtLayer, *logicHcal1AbsLayerEx;

  G4Material        *world_mat, *steel, *scint_mat, *alum, *copper;

  void DefineMaterials();
  G4VPhysicalVolume* ConstructDetector(); 
    
  void SetTiltViaNcross(const int ncross);

  PHG4HcalPrototype2DetectorMessenger*  fDetectorMessenger;

  double outer_tilt_angle;

  int active;
  int absorberactive;
  int layer;
  int blackhole;
  std::string detector_type;
  std::string superdetector;
  
};

#endif
