// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BLOCKSTEPPINGACTION_H
#define G4DETECTORS_PHG4BLOCKSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4VPhysicalVolume;
class PHG4BlockDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHParameters;
class PHG4Shower;

class PHG4BlockSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4BlockSteppingAction(PHG4BlockDetector *, const PHParameters *parameters);

  //! destructor
  virtual ~PHG4BlockSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step *, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode *);

 private:
  //! pointer to the detector
  PHG4BlockDetector *m_Detector;
  const PHParameters *m_Params;
  //! pointer to hit container
  PHG4HitContainer *m_HitContainer;
  PHG4Hit *m_Hit;
  PHG4Shower *m_SaveShower;
  G4VPhysicalVolume *m_SaveVolPre;
  G4VPhysicalVolume *m_SaveVolPost;
  int m_SaveTrackId;
  int m_SavePreStepStatus;
  int m_SavePostStepStatus;
  int m_ActiveFlag;
  int m_BlackHoleFlag;

  int m_UseG4StepsFlag;
};

#endif  //__G4PHPHYTHIAREADER_H__
