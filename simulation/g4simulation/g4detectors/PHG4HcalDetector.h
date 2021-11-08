// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4HCALDETECTOR_H
#define G4DETECTORS_PHG4HCALDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Region.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Types.hh>

#include <map>
#include <string>  // for string

class G4Material;
class G4LogicalVolume;
class G4UserSteppingAction;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;

class PHG4HcalDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4HcalDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const std::string& dnam, const int layer = 0);

  //! destructor
  ~PHG4HcalDetector(void) override
  {
  }

  //! construct
  void ConstructMe(G4LogicalVolume* world) override;

  void SetRadius(const G4double dbl) { radius = dbl * cm; }
  void SetLength(const G4double dbl) { length = dbl * cm; }
  void SetPosition(const G4double x, const G4double y, const G4double z)
  {
    xpos = x * cm;
    ypos = y * cm;
    zpos = z * cm;
  }

  void SetTilt(const double tilt) { _sciTilt = tilt; }
  void SetScintWidth(const double wid) { _sciWidth = wid * cm; }
  void SetNumScint(const int num) { _sciNum = num; }
  void SetScintPhi0(const G4double phi0) { _sciPhi0 = phi0; }  //  in units of sampling cells

  void SetThickness(const G4double dbl) { TrackerThickness = dbl * cm; }
  void SetMaterial(const std::string& mat) { material = mat; }
  void SetActive(const int i = 1) { active = i; }
  void SetAbsorberActive(const int i = 1) { absorberactive = i; }
  void SetDetectorType(const std::string& typ) { detector_type = typ; }
  int IsInCylinderActive(const G4VPhysicalVolume*);
  void SuperDetector(const std::string& name) { superdetector = name; }
  const std::string SuperDetector() const { return superdetector; }
  int get_Layer() const { return layer; }
  G4UserSteppingAction* GetSteppingAction() override
  {
    if (_region)
      return _region->GetRegionalSteppingAction();
    else
      return 0;
  }

  void Print(const std::string& what = "ALL") const override;

  static int INACTIVE;

 private:
  double GetLength(const double phi) const;
  G4Material* TrackerMaterial;
  G4double TrackerThickness;

  G4LogicalVolume* cylinder_logic;
  G4VPhysicalVolume* cylinder_physi;
  std::map<const G4VPhysicalVolume*, int> box_vol;

  G4double radius;
  G4double length;
  G4double xpos, ypos, zpos;
  G4double _sciTilt;
  G4double _sciWidth;
  G4int _sciNum;
  G4double _sciPhi0;  //  in units of sampling cells
  G4Region* _region;
  std::string material;
  int active;
  int absorberactive;
  int layer;
  std::string detector_type;
  std::string superdetector;
};

#endif
