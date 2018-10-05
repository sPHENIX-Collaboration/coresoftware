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

class PHG4Prototype2InnerHcalDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4Prototype2InnerHcalDetector(PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam);

  //! destructor
  virtual ~PHG4Prototype2InnerHcalDetector();

  //! construct
  virtual void Construct(G4LogicalVolume* world);

  virtual void Print(const std::string& what = "ALL") const;

  //!@name volume accessors
  //@{
  int IsInPrototype2InnerHcal(G4VPhysicalVolume*) const;
  //@}

  void SuperDetector(const std::string& name) { superdetector = name; }
  const std::string SuperDetector() const { return superdetector; }
  int get_Layer() const { return m_Layer; }

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

  int get_scinti_row_id(const std::string& volname);
  int get_steel_plate_id(const std::string& volname);

 protected:
  int ConstructInnerHcal(G4LogicalVolume* sandwich);
  int DisplayVolume(G4VSolid* volume, G4LogicalVolume* logvol, G4RotationMatrix* rotm = NULL);
  int DisplayVolume(G4LogicalVolume* volume, G4LogicalVolume* logvol, G4RotationMatrix* rotm = NULL);
  std::set<G4LogicalVolume*> m_ActiveVolumeSet;
  std::string superdetector;
  std::map<std::string, int> m_SteelPlateIdMap;
  std::map<std::string, int> m_ScintillatorIdMap;
  PHParameters* m_Params;
  G4LogicalVolume* m_InnerHcalSteelPlate;
  G4AssemblyVolume* m_InnerHcalAssembly;
  G4TwoVector m_SteelPlateCornerUpperLeft;
  G4TwoVector m_SteelPlateCornerUpperRight;
  G4TwoVector m_SteelPlateCornerLowerRight;
  G4TwoVector m_SteelPlateCornerLowerLeft;
  double m_ScintiUoneFrontSize;
  G4TwoVector m_ScintiUoneCornerUpperLeft;
  G4TwoVector m_ScintiUoneCornerUpperRight;
  G4TwoVector m_ScintiUoneCornerLowerRight;
  G4TwoVector m_ScintiUoneCornerLowerLeft;

  G4TwoVector m_ScintiU2CornerUpperLeft;
  G4TwoVector m_ScintiU2CornerUpperRight;
  G4TwoVector m_ScintiU2CornerLowerRight;
  G4TwoVector m_ScintiU2CornerLowerLeft;

  double m_ScintiT9DistanceToCorner;
  double m_ScintiT9FrontSize;
  G4TwoVector m_ScintiT9CornerUpperLeft;
  G4TwoVector m_ScintiT9CornerUpperRight;
  G4TwoVector m_ScintiT9CornerLowerRight;
  G4TwoVector m_ScintiT9CornerLowerLeft;

  double m_ScintiT10FrontSize;
  G4TwoVector m_ScintiT10CornerUpperLeft;
  G4TwoVector m_ScintiT10CornerUpperRight;
  G4TwoVector m_ScintiT10CornerLowerRight;
  G4TwoVector m_ScintiT10CornerLowerLeft;

  double m_ScintiT11FrontSize;
  G4TwoVector m_ScintiT11CornerUpperLeft;
  G4TwoVector m_ScintiT11CornerUpperRight;
  G4TwoVector m_ScintiT11CornerLowerRight;
  G4TwoVector m_ScintiT11CornerLowerLeft;

  double m_ScintiT12FrontSize;
  G4TwoVector m_ScintiT12CornerUpperLeft;
  G4TwoVector m_ScintiT12CornerUpperRight;
  G4TwoVector m_ScintiT12CornerLowerRight;
  G4TwoVector m_ScintiT12CornerLowerLeft;

  double m_ScintiX;
  double m_SteelZ;
  double m_SizeZ;
  double m_ScintiTileZ;
  double m_ScintiTileThickness;
  double m_ScintiBoxSmaller;
  double m_GapBetweenTiles;
  double m_ScintiGap;
  double m_DeltaPhi;
  double m_VolumeSteel;
  double m_VolumeScintillator;

  int m_NScintiPlates;
  int m_NSteelPlates;

  int m_ActiveFlag;
  int m_AbsorberActiveFlag;

  int m_Layer;
};

#endif
