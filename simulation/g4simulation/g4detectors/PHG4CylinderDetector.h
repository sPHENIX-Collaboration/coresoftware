#ifndef PHG4CylinderDetector_h
#define PHG4CylinderDetector_h

#include <g4main/PHG4Detector.h>

#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHG4CylinderSubsystem;
class PHParameters;

class PHG4CylinderDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4CylinderDetector(PHG4CylinderSubsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int layer = 0);

  //! destructor
  virtual ~PHG4CylinderDetector(void)
  {
  }

  //! construct
  void Construct(G4LogicalVolume *world);

  bool IsInCylinder(const G4VPhysicalVolume *) const;
  void SuperDetector(const std::string &name) { superdetector = name; }
  const std::string SuperDetector() const { return superdetector; }
  int get_Layer() const { return layer; }
 private:
  PHParameters *params;

  G4VPhysicalVolume *cylinder_physi;
  PHG4CylinderSubsystem *m_MySubSys;

  int layer;
  std::string superdetector;
};

#endif
