// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4SECTORDETECTOR_H
#define G4DETECTORS_PHG4SECTORDETECTOR_H

#include "PHG4SectorConstructor.h"

#include <g4main/PHG4Detector.h>

class G4LogicalVolume;
class PHG4SectorDisplayAction;
class PHG4SectorSubsystem;

class PHG4SectorDetector : public PHG4Detector, public PHG4Sector::PHG4SectorConstructor
{
 public:
  //! constructor
  PHG4SectorDetector(PHG4SectorSubsystem *subsys, PHCompositeNode *Node, const std::string &dnam = "SECTOR");

  //! destructor
  virtual ~PHG4SectorDetector(void)
  {
  }

  //! construct
  virtual void Construct(G4LogicalVolume *world);

  //!@name volume accessors
  //@{
  bool IsInSectorActive(G4VPhysicalVolume *);
  bool IsInSectorInactive(G4VPhysicalVolume *);
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
