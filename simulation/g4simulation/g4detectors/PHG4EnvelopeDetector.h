// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4ENVELOPEDETECTOR_H
#define G4DETECTORS_PHG4ENVELOPEDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4Types.hh>

#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;

class PHG4EnvelopeDetector : public PHG4Detector
{
 public:
  //Constructor
  PHG4EnvelopeDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, const std::string &dnam);

  //Destructor
  ~PHG4EnvelopeDetector() override {}

  //Construct
  void ConstructMe(G4LogicalVolume *world) override;

  //Volume accessors
  bool IsInEnvelope(G4VPhysicalVolume *) const;

  void SetPlace(G4double place_in_x, G4double place_in_y, G4double place_in_z)
  {
    _placeInX = place_in_x;
    _placeInY = place_in_y;
    _placeInZ = place_in_z;
  }

  void SetInnerRadius(G4double radius) { _innerRadius = radius; }
  void SetOuterRadius(G4double radius) { _outerRadius = radius; }
  void SetCylinderLength(G4double length) { _dZ_cyl = length; }

  void SetActive(const int i = 1) { _active = i; }
  int IsActive() const { return _active; }

  void SuperDetector(const std::string &name) { _superdetector = name; }
  const std::string SuperDetector() const { return _superdetector; }

  int get_Layer() const { return _layer; }

 private:
  G4double _placeInX;
  G4double _placeInY;
  G4double _placeInZ;
  G4double _innerRadius;
  G4double _outerRadius;
  G4double _dZ;
  G4double _dZ_cyl;
  G4double _sPhi;
  G4double _dPhi;

  G4String _materialCrystal;

  int _active;
  int _layer;

  std::string _superdetector;
};

#endif
