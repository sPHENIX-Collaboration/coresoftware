// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4JLEIC_G4JLEICVTXSTEPPINGACTION_H
#define G4JLEIC_G4JLEICVTXSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class G4JLeicVTXDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHParametersContainer;

class G4JLeicVTXSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  G4JLeicVTXSteppingAction(G4JLeicVTXDetector*, const PHParametersContainer*);

  //! destructor
  virtual ~G4JLeicVTXSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*);

 private:
  //! pointer to the detector
  G4JLeicVTXDetector* m_Detector;
  //! pointer to hit container
  PHG4HitContainer* m_HitContainer;
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

#endif  // G4JLEIC_G4JLEICVTXSTEPPINGACTION_H
