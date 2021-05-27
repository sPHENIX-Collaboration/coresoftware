// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CONEDETECTOR_H
#define G4DETECTORS_PHG4CONEDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <string>  // for string

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4ConeDisplayAction;
class PHG4Subsystem;
class PHParameters;

class PHG4ConeDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4ConeDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam, const int lyr = 0);

  //! destructor
  ~PHG4ConeDetector(void) override
  {
  }

  //! construct
  void ConstructMe(G4LogicalVolume* world) override;

  //!@name volume accessors
  //@{
  bool IsInConeActive(G4VPhysicalVolume*);
  //@}

  void SuperDetector(const std::string& name) { superdetector = name; }
  const std::string SuperDetector() const { return superdetector; }
  int get_Layer() const { return layer; }

 private:
  PHParameters* m_Params = nullptr;

  G4VPhysicalVolume* m_ConePhysVol = nullptr;
  PHG4ConeDisplayAction* m_DisplayAction = nullptr;

  int layer = -9999;
  std::string superdetector;
};

#endif
