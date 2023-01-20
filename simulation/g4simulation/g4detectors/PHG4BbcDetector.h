// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BBCDETECTOR_H
#define G4DETECTORS_PHG4BBCDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <cmath>
#include <set>
#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHG4BbcDisplayAction;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;

class PHG4BbcDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4BbcDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *params, const std::string &dnam = "BBC");

  //! destructor
  ~PHG4BbcDetector() override
  {
  }

  //! construct BBC/MBD
  void ConstructMe(G4LogicalVolume *world) override;
  void ConstructSupport(G4LogicalVolume *world);

  void Print(const std::string &what = "ALL") const override;

  int IsInBbc(G4VPhysicalVolume *) const;

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }

  PHG4BbcDisplayAction *GetDisplayAction() { return m_DisplayAction; }

 protected:
  PHG4BbcDisplayAction *m_DisplayAction = nullptr;

  PHParameters *m_Params = nullptr;

  int m_ActiveFlag = 0;
  int m_SupportActiveFlag = 0;

  float m_bbcz = NAN;  // z-location of mid-point of quartz ckov crystals

  std::set<G4LogicalVolume *> m_PhysLogicalVolSet;
  std::set<G4LogicalVolume *> m_SupportLogicalVolSet;

  std::string m_SuperDetector;
};

#endif
