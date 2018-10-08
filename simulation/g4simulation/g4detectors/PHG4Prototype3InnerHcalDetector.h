// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4DETECTORS_PHG4PROTOTYPE3INNERHCALDETECTOR_H
#define G4DETECTORS_PHG4PROTOTYPE3INNERHCALDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4TwoVector.hh>

#include <map>
#include <set>
#include <string>

class G4AssemblyVolume;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;
class PHParameters;

class PHG4Prototype3InnerHcalDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4Prototype3InnerHcalDetector(PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  virtual ~PHG4Prototype3InnerHcalDetector();
  //! construct method called by G4
  void Construct(G4LogicalVolume *world);

  //! print detector info
  void Print(const std::string &what = "ALL") const;

  //!@name volume accessors
  //@{
  int IsInPrototype3InnerHcal(G4VPhysicalVolume *) const;
  //@}

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
  int get_Layer() const { return m_Layer; }

  int get_scinti_row_id(const std::string &volname);
  int get_steel_plate_id(const std::string &volname);

 private:
  G4LogicalVolume *ConstructSteelPlate(G4LogicalVolume *hcalenvelope);
  G4LogicalVolume *ConstructScintillatorBoxHiEta(G4LogicalVolume *hcalenvelope);
  G4LogicalVolume *ConstructScintiTile9(G4LogicalVolume *hcalenvelope);
  G4LogicalVolume *ConstructScintiTile10(G4LogicalVolume *hcalenvelope);
  G4LogicalVolume *ConstructScintiTile11(G4LogicalVolume *hcalenvelope);
  G4LogicalVolume *ConstructScintiTile12(G4LogicalVolume *hcalenvelope);
  double GetScintiAngle();
  int ConstructInnerHcal(G4LogicalVolume *sandwich);

  PHParameters *m_params;
  G4LogicalVolume *m_InnerHcalSteelPlate;
  G4AssemblyVolume *m_InnerHcalAssembly;
  G4LogicalVolume *m_scintibox;
  std::set<G4LogicalVolume *> m_ActiveVolumeSet;
  G4TwoVector m_SteelPlateCornerUpperLeft;
  G4TwoVector m_SteelPlateCornerUpperRight;
  G4TwoVector m_SteelPlateCornerLowerRight;
  G4TwoVector m_SteelPlateCornerLowerLeft;

  double m_ScintiTile9DistanceToCorner;
  double m_ScintiTile9FrontSize;
  G4TwoVector m_ScintiTile9CornerUpperLeft;
  G4TwoVector m_ScintiTile9CornerUpperRight;
  G4TwoVector m_ScintiTile9CornerLowerRight;
  G4TwoVector m_ScintiTile9CornerLowerLeft;

  double m_ScintiTile10FrontSize;
  G4TwoVector m_ScintiTile10CornerUpperLeft;
  G4TwoVector m_ScintiTile10CornerUpperRight;
  G4TwoVector m_ScintiTile10CornerLowerRight;
  G4TwoVector m_ScintiTile10CornerLowerLeft;

  double m_ScintiTile11FrontSize;
  G4TwoVector m_ScintiTile11CornerUpperLeft;
  G4TwoVector m_ScintiTile11CornerUpperRight;
  G4TwoVector m_ScintiTile11CornerLowerRight;
  G4TwoVector m_ScintiTile11CornerLowerLeft;

  double m_ScintiTile12FrontSize;
  G4TwoVector m_ScintiTile12CornerUpperLeft;
  G4TwoVector m_ScintiTile12CornerUpperRight;
  G4TwoVector m_ScintiTile12CornerLowerRight;
  G4TwoVector m_ScintiTile12CornerLowerLeft;

  double m_ScintiX;
  double m_SteelZ;
  double m_ScintiTileZ;
  double m_ScintiTileThickness;
  double m_GapBetweenTiles;
  double m_ScintiGap;
  double m_DeltaPhi;
  double m_VolumeSteel;
  double m_VolumeScintillator;

  G4TwoVector m_ScintiCornerLowerLeft;
  G4TwoVector m_ScintiCornerLowerRight;
  // leave this in in case we ever need those coordinates
  // G4TwoVector m_ScintiCornerUpperRight;
  // G4TwoVector m_ScintiCornerUpperLeft;

  int m_NumScintiPlates;
  int m_NumSteelPlates;

  int m_Active;
  int m_AbsorberActive;

  int m_Layer;
  std::string m_SuperDetector;
  std::map<std::string,int> m_SteelPlateIdMap;
  std::map<std::string,int> m_ScintillatorIdMap;
};

#endif  // G4DETECTORS_PHG4PROTOTYPE3INNERHCALDETECTOR_H
