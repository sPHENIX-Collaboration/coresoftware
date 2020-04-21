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

class BeamLineMagnetDetector : public PHG4Detector
{
 public:
  //! constructor
  BeamLineMagnetDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int layer = 0);

  //! destructor
  virtual ~BeamLineMagnetDetector(void)
  {
  }

  //! construct
  void ConstructMe(G4LogicalVolume *world);

  int IsInBeamLineMagnet(const G4VPhysicalVolume *) const;
  void SuperDetector(const std::string &name) { superdetector = name; }
  const std::string SuperDetector() const { return superdetector; }
  int get_Layer() const { return layer; }

 private:
  PHParameters *params;

  G4VPhysicalVolume *magnet_physi;
  G4VPhysicalVolume *magnet_iron_physi;
  BeamLineMagnetDisplayAction *m_DisplayAction;

  int layer;
  std::string superdetector;
};

#endif  //  G4DETECTORS_BEAMLINEMAGNETMAGNETDETECTOR_H
