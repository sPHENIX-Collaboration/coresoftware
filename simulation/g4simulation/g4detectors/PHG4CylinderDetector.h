// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERDETECTOR_H
#define G4DETECTORS_PHG4CYLINDERDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHG4CylinderDisplayAction;
class PHG4CylinderSubsystem;
class PHParameters;

class PHG4CylinderDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4CylinderDetector(PHG4CylinderSubsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int layer = 0);

  //! destructor
  virtual ~PHG4CylinderDetector()
  {
  }

  //! construct
  void Construct(G4LogicalVolume *world);

  bool IsInCylinder(const G4VPhysicalVolume *) const;
  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
  int get_Layer() const { return m_Layer; }

 private:
  PHParameters *m_Params;

  G4VPhysicalVolume *m_CylinderPhysicalVolume;
  PHG4CylinderDisplayAction *m_DisplayAction;

  int m_Layer;
  std::string m_SuperDetector;
};

#endif
