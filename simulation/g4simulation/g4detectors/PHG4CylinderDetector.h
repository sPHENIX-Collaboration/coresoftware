// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERDETECTOR_H
#define G4DETECTORS_PHG4CYLINDERDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4CylinderDisplayAction;
class PHG4Subsystem;
class PHParameters;

class PHG4CylinderDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4CylinderDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int layer = 0);

  //! destructor
  ~PHG4CylinderDetector() override
  {
  }

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

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
