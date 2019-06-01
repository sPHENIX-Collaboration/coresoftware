// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4DETECTORS_PHG4INNERHCALSTEPPINGACTION_H
#define G4DETECTORS_PHG4INNERHCALSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4InnerHcalDetector;
class PHParameters;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;

class PHG4InnerHcalSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4InnerHcalSteppingAction(PHG4InnerHcalDetector *, const PHParameters *parameters);

  //! destructor
  virtual ~PHG4InnerHcalSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step *, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode *);

 private:
  //! pointer to the detector
  PHG4InnerHcalDetector *m_Detector;

  //! pointer to hit container
  PHG4HitContainer *m_Hits;
  PHG4HitContainer *m_AbsorberHits;
  PHG4Hit *m_Hit;
  const PHParameters *m_Params;
  PHG4HitContainer *m_SaveHitContainer;
  PHG4Shower *m_SaveShower;
  G4VPhysicalVolume *m_SaveVolPre;
  G4VPhysicalVolume *m_SaveVolPost;
  int m_SaveTrackId;
  int m_SavePreStepStatus;
  int m_SavePostStepStatus;
  // since getting parameters is a map search we do not want to
  // do this in every step, the parameters used are cached
  // in the following variables
  int m_IsActive;
  int m_IsBlackHole;
  int m_LightScintModel;
};

#endif  // G4DETECTORS_PHG4INNERHCALSTEPPINGACTION_H
