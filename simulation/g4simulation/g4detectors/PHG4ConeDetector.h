// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CONEDETECTOR_H
#define G4DETECTORS_PHG4CONEDETECTOR_H

#include "g4main/PHG4Detector.h"

#include <Geant4/G4Region.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Types.hh>

#include <string>  // for string

class G4Cons;
class G4LogicalVolume;
class G4Material;
class G4UserSteppingAction;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;

class PHG4ConeDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4ConeDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const std::string& dnam, const int lyr = 0);

  //! destructor
  virtual ~PHG4ConeDetector(void)
  {
  }

  //! construct
  virtual void ConstructMe(G4LogicalVolume* world);

  //!@name volume accessors
  //@{
  bool IsInConeActive(G4VPhysicalVolume*);
  bool IsInConeInactive(G4VPhysicalVolume*);
  //@}

  //!set inner and outter radius1
  void SetR1(const G4double min, const G4double max)
  {
    rMin1 = min * cm;
    rMax1 = max * cm;
  }

  //!set inner and outter radius2
  void SetR2(const G4double min, const G4double max)
  {
    rMin2 = min * cm;
    rMax2 = max * cm;
  }

  //! set length in Z
  void SetZlength(const G4double a)
  {
    dZ = a * cm;
  }

  //! set phi offset and extention
  void SetPhi(const G4double a, const G4double b)
  {
    sPhi = a;
    dPhi = g;
  }

  void SetMaterial(const std::string& mat) { material = mat; }
  void SetPlaceZ(const G4double place_z) { place_in_z = place_z * cm; }
  void SetPlace(const G4double place_x, const G4double place_y, const G4double place_z)
  {
    place_in_x = place_x * cm;
    place_in_y = place_y * cm;
    place_in_z = place_z * cm;
  }
  void SetZRot(const G4double z_angle) { z_rot = z_angle * rad; }
  virtual G4UserSteppingAction* GetSteppingAction()
  {
    if (_region)
      return _region->GetRegionalSteppingAction();
    else
      return 0;
  }
  void SetActive(const int i = 1) { active = i; }
  void SuperDetector(const std::string& name) { superdetector = name; }
  const std::string SuperDetector() const { return superdetector; }
  int get_Layer() const { return layer; }

 private:
  G4Material* TrackerMaterial;

  G4Cons* block_solid;
  G4LogicalVolume* block_logic;
  G4VPhysicalVolume* block_physi;

  G4String material;
  G4double place_in_x;
  G4double place_in_y;
  G4double place_in_z;
  G4double rMin1;
  G4double rMax1;
  G4double rMin2;
  G4double rMax2;
  G4double dZ;
  G4double sPhi;
  G4double dPhi;
  G4double z_rot;

  G4Region* _region;
  int active;
  int layer;
  std::string superdetector;
};

#endif
