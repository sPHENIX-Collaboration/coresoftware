#ifndef PHG4HcalTestBeamDetector_h
#define PHG4HcalTestBeamDetector_h

#include "g4main/PHG4Detector.h"

#include <Geant4/globals.hh>
#include <Geant4/G4Types.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4RotationMatrix.hh>

#include <map>
#include <vector>

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;

class HcalGeom
{

 public:
  HcalGeom(const std::string &nam);
  virtual ~HcalGeom() {}
  void InitInner();
  void InitOuter();
  void CalculateGeometryInner();
  void CalculateGeometryOuter();
  void CalculateGeometry();

  std::string name;
  int num_sandwiches;
  double steel_plate_x;
  double steel_plate_y_in;
  double steel_plate_y_out;
  double steel_plate_z;
  double scintillator_x;
  double scintillator_y;
  double scintillator_z;
  double scintillator_y_vertical_projection;
  double tilt_angle;
  double inner_radius;
  double scinti_gap;

 protected:
  HcalGeom() {}
};


class PHG4HcalTestBeamDetector: public PHG4Detector
{

  public:

  //! constructor
  PHG4HcalTestBeamDetector( PHCompositeNode *Node, const std::string &dnam="BLOCK", const int lyr = 0 );

  //! destructor
  virtual ~PHG4HcalTestBeamDetector();

  //! construct
  virtual void Construct( G4LogicalVolume* world );

  //!@name volume accessors
  //@{
  int IsInHcalTestBeam(G4VPhysicalVolume*) const;
  //@}

  void SetPlaceZ(const G4double place_z) {place_in_z = place_z*cm;}
  void SetPlace(const G4double place_x, const G4double place_y, const G4double place_z)
  {
    place_in_x = place_x*cm;
    place_in_y = place_y*cm;
    place_in_z = place_z*cm;
  }
  void SetXRot(const G4double angle) {x_rot = angle*rad;}
  void SetYRot(const G4double angle) {y_rot = angle*rad;}
  void SetZRot(const G4double angle) {z_rot = angle*rad;}
  void SetActive(const int i = 1) {active = i;}
  void SetAbsorberActive(const int i = 1) {absorberactive = i;}
  int IsActive() const {return active;}
  void SuperDetector(const std::string &name) {superdetector = name;}
  const std::string SuperDetector() const {return superdetector;}
  int get_Layer() const {return layer;}

  void BlackHole(const int i=1) {blackhole = i;}
  int IsBlackHole() const {return blackhole;}

  private:

  void CalculateGeometry();
  int ConstructOuterHcalTower(G4LogicalVolume* tower);
  int ConstructInnerHcalTower(G4LogicalVolume* tower);
  int ConstructOuterHcal(G4LogicalVolume* tower);
  int ConstructInnerHcal(G4LogicalVolume* tower);
  int ConstructOuterSandwichVolume(G4LogicalVolume* sandwich);
  int ConstructInnerSandwichVolume(G4LogicalVolume* sandwich);
  int ConstructInnerAbsorber(G4LogicalVolume* sandwich);
  int ConstructAbsorber(G4LogicalVolume* sandwich);
  int ConstructTest(G4LogicalVolume* sandwich);

  int DisplayVolume(G4VSolid *volume,  G4LogicalVolume* logvol, G4RotationMatrix* rotm=NULL);

  G4VPhysicalVolume *inner_scinti;
  G4VPhysicalVolume *inner_absorber;
  G4VPhysicalVolume *outer_scinti;
  G4VPhysicalVolume *outer_absorber;
  double outer_plate_x;
  double outer_plate_z;
  G4double outer_steel_x;
  G4double outer_steel_y_in; // inner diameter
  G4double outer_steel_y_out; // outer diameter
  G4double outer_steel_z;
  G4double inner_steel_x;
  G4double inner_steel_y_in; // inner diameter
  G4double inner_steel_y_out; // outer diameter
  G4double inner_steel_z;
  G4double outer_sc_dimension[3];
  G4double inner_sc_dimension[3];
  G4double place_in_x;
  G4double place_in_y;
  G4double place_in_z;
  G4double x_rot;
  G4double y_rot;
  G4double z_rot;
  double outer_tilt_angle;
  double inner_tilt_angle;
  double inner_radius;
  double outer_radius;
  int num_sandwiches;
  double inner_hcal_inner_radius;
  double outer_hcal_inner_radius;
  int active;
  int absorberactive;
  int layer;
  int blackhole;
  std::string detector_type;
  std::string superdetector;
  HcalGeom *OuterGeom;
  HcalGeom *InnerGeom;
};

#endif
