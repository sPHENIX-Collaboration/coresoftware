// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_BEAMLINEMAGNETMAGNETDETECTOR_H
#define G4DETECTORS_BEAMLINEMAGNETMAGNETDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <string>

class BeamLineMagnetDisplayAction;
class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;
class G4MagneticField;

class BeamLineMagnetDetector : public PHG4Detector
{
 public:
  //! constructor
  BeamLineMagnetDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int magnetid = 0);

  //! destructor
  ~BeamLineMagnetDetector(void) override
  {
  }

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  //! Optional PostConstruction call after all geometry is constructed
  void PostConstruction() override;

  int IsInBeamLineMagnet(const G4VPhysicalVolume *) const;
  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
  int get_MagnetId() const { return m_MagnetId; }

 private:
  PHParameters *m_Params = nullptr;

  G4VPhysicalVolume *magnet_physi = nullptr;
  G4VPhysicalVolume *magnet_iron_physi = nullptr;
  G4VPhysicalVolume *magnet_core_physi = nullptr;
  BeamLineMagnetDisplayAction *m_DisplayAction = nullptr;
  G4LogicalVolume *m_magnetFieldLogic = nullptr;
  G4MagneticField *m_magField = nullptr;

  int m_MagnetId = -1;

  std::string m_SuperDetector = "NONE";
};

#endif  //  G4DETECTORS_BEAMLINEMAGNETMAGNETDETECTOR_H
