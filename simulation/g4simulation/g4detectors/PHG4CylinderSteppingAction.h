// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERSTEPPINGACTION_H
#define G4DETECTORS_PHG4CYLINDERSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <string>

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4CylinderDetector;
class PHG4CylinderSubsystem;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;
class PHParameters;

class PHG4CylinderSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4CylinderSteppingAction(PHG4CylinderSubsystem *subsys, PHG4CylinderDetector *detector, const PHParameters *parameters);

  //! destructor
  ~PHG4CylinderSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step *, bool) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode *) override;

  void SaveLightYield(const int i = 1) { m_SaveLightYieldFlag = i; }

  // needed for hit position crosschecks, if this volume is inside
  // another volume the absolut hit coordinates in our G4Hits and
  // the local coordinates differ, so checking against our place in z
  // goes wrong
  bool hasMotherSubsystem() const;

  // this is just needed for use as reference plane for projections
  // this is the only detector using this - there is no need to add
  // this to our parameters
  void SaveAllHits(bool i = true) { m_SaveAllHitsFlag = i; }

  void HitNodeName(const std::string &name) { m_HitNodeName = name; }

 private:
  //! pointer to the Subsystem
  PHG4CylinderSubsystem *m_Subsystem;
  //! pointer to the detector
  PHG4CylinderDetector *m_Detector;
  const PHParameters *m_Params;

  //! pointer to hit container
  PHG4HitContainer *m_HitContainer;
  PHG4Hit *m_Hit;
  PHG4Shower *m_SaveShower;
  G4VPhysicalVolume *m_SaveVolPre;
  G4VPhysicalVolume *m_SaveVolPost;
  bool m_SaveAllHitsFlag = false;
  int m_SaveLightYieldFlag;
  int m_SaveTrackId;
  int m_SavePreStepStatus;
  int m_SavePostStepStatus;
  int m_ActiveFlag;
  int m_BlackHoleFlag;
  int m_UseG4StepsFlag;
  double m_Zmin;
  double m_Zmax;
  double m_Tmin;
  double m_Tmax;
  std::string m_HitNodeName;
};

#endif
