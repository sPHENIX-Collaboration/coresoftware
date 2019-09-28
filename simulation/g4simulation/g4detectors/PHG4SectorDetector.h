// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4SECTORDETECTOR_H
#define G4DETECTORS_PHG4SECTORDETECTOR_H

#include "PHG4SectorConstructor.h"

#include <g4main/PHG4Detector.h>

#include <string>  // for string

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4SectorDisplayAction;
class PHG4Subsystem;

class PHG4SectorDetector : public PHG4Detector, public PHG4Sector::PHG4SectorConstructor
{
 public:
  //! constructor
  PHG4SectorDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, const std::string &dnam);

  //! destructor
  virtual ~PHG4SectorDetector(void)
  {
  }

  //! construct
  virtual void ConstructMe(G4LogicalVolume *world);

  //!@name volume accessors
  //@{
  bool IsInSectorActive(G4VPhysicalVolume *physvol);
  //@}

  void SuperDetector(const std::string &name) { superdetector = name; }
  const std::string SuperDetector() const { return superdetector; }

  virtual void OverlapCheck(const bool chk = true)
  {
    PHG4Detector::OverlapCheck(chk);
    PHG4SectorConstructor::OverlapCheck(chk);
  }

 private:
  PHG4SectorDisplayAction *m_DisplayAction;

  std::string superdetector;
};

#endif
