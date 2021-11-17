// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BLOCKDETECTOR_H
#define G4DETECTORS_PHG4BLOCKDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <string>  // for string

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4BlockDisplayAction;
class PHG4Subsystem;
class PHParameters;

class PHG4BlockDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4BlockDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int lyr = 0);

  //! destructor
  ~PHG4BlockDetector(void) override
  {
  }

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  //!@name volume accessors
  //@{
  bool IsInBlock(G4VPhysicalVolume *) const;
  //@}

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
  int get_Layer() const { return m_Layer; }

 private:
  PHParameters *m_Params;

  G4VPhysicalVolume *m_BlockPhysi;
  PHG4BlockDisplayAction *m_DisplayAction;

  int m_Layer;
  std::string m_SuperDetector;
};

#endif
