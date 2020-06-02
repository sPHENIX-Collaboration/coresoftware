// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHG4MICROMEGASSTEPPINGACTION_H
#define PHG4MICROMEGASSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class PHG4MicromegasDetector;

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Hit;
class PHG4HitContainer;
class PHParameters;

class PHG4MicromegasSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4MicromegasSteppingAction(PHG4MicromegasDetector*, const PHParameters* parameters);

  //! destructor
  virtual ~PHG4MicromegasSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*);

 private:
  //! pointer to the detector
  PHG4MicromegasDetector* m_Detector;
  const PHParameters* m_Params;
  //! pointer to hit container
  PHG4HitContainer* m_HitContainer;
  PHG4Hit* m_Hit;
  PHG4HitContainer* m_SaveHitContainer;
  G4VPhysicalVolume* m_SaveVolPre;
  G4VPhysicalVolume* m_SaveVolPost;

  int m_SaveTrackId;
  int m_SavePreStepStatus;
  int m_SavePostStepStatus;
  int m_ActiveFlag;
  int m_BlackHoleFlag;
  double m_EdepSum;
  double m_EionSum;
};

#endif // MICROMEGASSTEPPINGACTION_H
