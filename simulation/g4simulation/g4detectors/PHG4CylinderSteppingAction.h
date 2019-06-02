// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERSTEPPINGACTION_H
#define G4DETECTORS_PHG4CYLINDERSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4CylinderDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;
class PHParameters;

class PHG4CylinderSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4CylinderSteppingAction(PHG4CylinderDetector *, const PHParameters *parameters);

  //! destructor
  virtual ~PHG4CylinderSteppingAction();

  //! stepping action
  bool UserSteppingAction(const G4Step *, bool);

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode *);

  void SaveLightYield(const int i = 1) { m_SaveLightYieldFlag = i; }

 private:
  //! pointer to the detector
  PHG4CylinderDetector *m_Detector;

  const PHParameters *m_Params;

  //! pointer to hit container
  PHG4HitContainer *m_HitContainer;
  PHG4Hit *m_Hit;
  PHG4Shower *m_SaveShower;
  G4VPhysicalVolume *m_SaveVolPre;
  G4VPhysicalVolume *m_SaveVolPost;
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
};

#endif
