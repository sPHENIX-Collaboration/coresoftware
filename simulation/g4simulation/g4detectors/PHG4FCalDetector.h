// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4FCALDETECTOR_H
#define G4DETECTORS_PHG4FCALDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Region.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4Types.hh>

#include <map>

class G4Material;
class G4Box;
class G4LogicalVolume;
class G4UserSteppingAction;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;

class PHG4FCalDetector : public PHG4Detector
{
 public:
  PHG4FCalDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const std::string& name);

  virtual ~PHG4FCalDetector(void) {}

  virtual void ConstructMe(G4LogicalVolume* world);

  virtual G4UserSteppingAction* GetSteppingAction()
  {
    if (_region)
      return _region->GetRegionalSteppingAction();
    else
      return nullptr;
  }

  bool isInScintillator(G4VPhysicalVolume* volume);
  int getScintillatorLayer(G4VPhysicalVolume* volume);
  unsigned int computeIndex(unsigned int layer, G4double x, G4double y, G4double z, G4double& xcenter, G4double& ycenter, G4double& zcenter);

 private:
  G4double length, absorber_thickness, scintillator_thickness;
  unsigned int nlayers, segments_per_column, segments_per_thickness;
  G4double z_position;
  G4double layer_separation;

  G4Material* SetMaterial(G4String);

  G4Material* AbsorberMaterial;
  G4Material* ScintillatorMaterial;

  std::map<unsigned int, G4Box*> absorber_solid_;
  std::map<unsigned int, G4LogicalVolume*> absorber_logic_;
  std::map<unsigned int, G4VPhysicalVolume*> absorber_physi_;

  std::map<unsigned int, G4Box*> scintillator_solid_;
  std::map<unsigned int, G4LogicalVolume*> scintillator_logic_;
  std::map<unsigned int, G4VPhysicalVolume*> scintillator_physi_;

  G4Region* _region;
};

#endif
