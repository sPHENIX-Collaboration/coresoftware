// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BBCSTEPPINGACTION_H
#define G4DETECTORS_PHG4BBCSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class PHG4BbcDetector;

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Hit;
class PHG4HitContainer;
class PHParameters;

class PHG4BbcSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4BbcSteppingAction(PHG4BbcDetector*, const PHParameters*);

  //! destructor
  virtual ~PHG4BbcSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*);

 private:
  //! pointer to the detector
  PHG4BbcDetector* m_Detector;

  //! pointer to hit container
  PHG4HitContainer* m_HitContainer;
  PHG4Hit* m_Hit;
  PHG4HitContainer *m_SaveHitContainer;

  G4VPhysicalVolume *m_SaveVolPre;
  G4VPhysicalVolume *m_SaveVolPost;

  int m_SaveTrackId;
  int m_SavePreStepStatus;
  int m_SavePostStepStatus;
  //int m_ActiveFlag;
  //int m_BlackHoleFlag;
  double m_EdepSum;
  double m_EionSum;

  //! pointer to parameters
  //const PHParameters* m_Params;

};

#endif  // PHG4BbcSteppingAction_h__
