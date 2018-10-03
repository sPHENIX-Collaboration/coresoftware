// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4DETECTORS_PHG4PROTOTYPE2INNERHCALDETECTOR_H
#define G4DETECTORS_PHG4PROTOTYPE2INNERHCALDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <Geant4/G4TwoVector.hh>

#include <map>
#include <set>

class G4AssemblyVolume;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;
class PHParameters;

class PHG4Prototype2InnerHcalDetector: public PHG4Detector
{

  public:

  //! constructor
 PHG4Prototype2InnerHcalDetector( PHCompositeNode *Node,  PHParameters *parameters, const std::string &dnam);

  //! destructor
 virtual ~PHG4Prototype2InnerHcalDetector();

  //! construct
  virtual void Construct( G4LogicalVolume* world );

  virtual void Print(const std::string &what = "ALL") const;

  //!@name volume accessors
  //@{
  int IsInPrototype2InnerHcal(G4VPhysicalVolume*) const;
  //@}

  void SuperDetector(const std::string &name) {superdetector = name;}
  const std::string SuperDetector() const {return superdetector;}
  int get_Layer() const {return layer;}

  G4LogicalVolume* ConstructSteelPlate(G4LogicalVolume* hcalenvelope);
  G4LogicalVolume* ConstructScintillatorBox(G4LogicalVolume* hcalenvelope);
  G4LogicalVolume* ConstructScintillatorBoxHiEta(G4LogicalVolume* hcalenvelope);
  G4LogicalVolume* ConstructScintiTileU1(G4LogicalVolume* hcalenvelope);
  G4LogicalVolume* ConstructScintiTileU2(G4LogicalVolume* hcalenvelope);
  G4LogicalVolume* ConstructScintiTile9(G4LogicalVolume* hcalenvelope);
  G4LogicalVolume* ConstructScintiTile10(G4LogicalVolume* hcalenvelope);
  G4LogicalVolume* ConstructScintiTile11(G4LogicalVolume* hcalenvelope);
  G4LogicalVolume* ConstructScintiTile12(G4LogicalVolume* hcalenvelope);
  double GetScintiAngle();

  int get_scinti_row_id(const std::string &volname);
  int get_steel_plate_id(const std::string &volname);

  protected:
  int ConstructInnerHcal(G4LogicalVolume* sandwich);
  int DisplayVolume(G4VSolid *volume,  G4LogicalVolume* logvol, G4RotationMatrix* rotm=NULL);
  int DisplayVolume(G4LogicalVolume *volume,  G4LogicalVolume* logvol, G4RotationMatrix* rotm=NULL);
  PHParameters *params;
  G4LogicalVolume *m_InnerHcalSteelPlate;
  G4AssemblyVolume *m_InnerHcalAssembly;
  std::set<G4LogicalVolume *> m_ActiveVolumeSet;
  G4TwoVector steel_plate_corner_upper_left;
  G4TwoVector steel_plate_corner_upper_right;
  G4TwoVector steel_plate_corner_lower_right;
  G4TwoVector steel_plate_corner_lower_left;
  double scinti_u1_front_size;
  G4TwoVector scinti_u1_corner_upper_left;
  G4TwoVector scinti_u1_corner_upper_right;
  G4TwoVector scinti_u1_corner_lower_right;
  G4TwoVector scinti_u1_corner_lower_left;

  G4TwoVector scinti_u2_corner_upper_left;
  G4TwoVector scinti_u2_corner_upper_right;
  G4TwoVector scinti_u2_corner_lower_right;
  G4TwoVector scinti_u2_corner_lower_left;

  double scinti_t9_distance_to_corner;
  double scinti_t9_front_size;
  G4TwoVector scinti_t9_corner_upper_left;
  G4TwoVector scinti_t9_corner_upper_right;
  G4TwoVector scinti_t9_corner_lower_right;
  G4TwoVector scinti_t9_corner_lower_left;

  double scinti_t10_front_size;
  G4TwoVector scinti_t10_corner_upper_left;
  G4TwoVector scinti_t10_corner_upper_right;
  G4TwoVector scinti_t10_corner_lower_right;
  G4TwoVector scinti_t10_corner_lower_left;

  double scinti_t11_front_size;
  G4TwoVector scinti_t11_corner_upper_left;
  G4TwoVector scinti_t11_corner_upper_right;
  G4TwoVector scinti_t11_corner_lower_right;
  G4TwoVector scinti_t11_corner_lower_left;

  double scinti_t12_front_size;
  G4TwoVector scinti_t12_corner_upper_left;
  G4TwoVector scinti_t12_corner_upper_right;
  G4TwoVector scinti_t12_corner_lower_right;
  G4TwoVector scinti_t12_corner_lower_left;

  double scinti_x;
  double steel_z;
  double size_z;
  double scinti_tile_z;
  double scinti_tile_thickness;
  double scinti_box_smaller;
  double gap_between_tiles;
  double scinti_gap;
  double deltaphi;
  double volume_steel;
  double volume_scintillator;

  int n_scinti_plates;
  int n_steel_plates;

  int m_ActiveFlag;
  int m_AbsorberActiveFlag;

  int layer;
//  std::string detector_type;
  std::string superdetector;
  std::map<std::string,int> m_SteelPlateIdMap;
  std::map<std::string,int> m_ScintillatorIdMap;

};

#endif
