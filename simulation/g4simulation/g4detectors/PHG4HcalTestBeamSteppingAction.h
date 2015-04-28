#ifndef PHG4VHcalTestBeamSteppingAction_h
#define PHG4VHcalTestBeamSteppingAction_h

#include "g4main/PHG4SteppingAction.h"

class PHG4HcalTestBeamDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4HcalTestBeamSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4HcalTestBeamSteppingAction( PHG4HcalTestBeamDetector* );

  //! destroctor
  virtual ~PHG4HcalTestBeamSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  private:

  //! pointer to the detector
  PHG4HcalTestBeamDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4HitContainer * absorberhits_;
  PHG4Hit *hit;
};


#endif // PHG4HcalTestBeamSteppingAction_h
