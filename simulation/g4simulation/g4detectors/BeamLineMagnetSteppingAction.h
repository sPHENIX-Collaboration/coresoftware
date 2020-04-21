// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BEAMLINEMAGNETSTEPPINGACTION_H
#define BEAMLINEMAGNETSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class BeamLineMagnetDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHParameters;

class BeamLineMagnetSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  BeamLineMagnetSteppingAction(BeamLineMagnetDetector*, const PHParameters *parameters);

  //! destructor
  virtual ~BeamLineMagnetSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*);

 private:
  //! pointer to the detector
  BeamLineMagnetDetector* m_Detector;
  //! pointer to hit container
  PHG4HitContainer* m_HitContainer;
  PHG4HitContainer* m_AbsorberHitContainer;
  PHG4Hit* m_Hit;
  PHG4HitContainer* m_SaveHitContainer;

  G4VPhysicalVolume* m_SaveVolPre;
  G4VPhysicalVolume* m_SaveVolPost;
  int m_SaveTrackId;
  int m_SavePreStepStatus;
  int m_SavePostStepStatus;
  double m_EdepSum;
  double m_EionSum;
};

#endif  // BEAMLINEMAGNETSTEPPINGACTION_H
