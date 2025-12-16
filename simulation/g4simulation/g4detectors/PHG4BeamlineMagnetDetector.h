// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BEAMLINEMAGNETDETECTOR_H
#define G4DETECTORS_PHG4BEAMLINEMAGNETDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;

class PHG4BeamlineMagnetDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4BeamlineMagnetDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int lyr = 0);

  //! destructor
  ~PHG4BeamlineMagnetDetector() override = default;

  //! construct
  void ConstructMe(G4LogicalVolume *logicMother) override;

  bool IsInBeamlineMagnet(const G4VPhysicalVolume *) const;
  void SuperDetector(const std::string &name) { superdetector = name; }
  const std::string &SuperDetector() const { return superdetector; }
  int get_Layer() const { return layer; }

 private:
  PHParameters *params {nullptr};

  G4VPhysicalVolume *magnet_physi {nullptr};
  G4VPhysicalVolume *cylinder_physi {nullptr};

  int layer {0};
  std::string superdetector;
};

#endif
