#ifndef PHG4VCEmcTestBeamSteppingAction_h
#define PHG4VCEmcTestBeamSteppingAction_h

#include "g4main/PHG4SteppingAction.h"

class PHG4CEmcTestBeamDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4CEmcTestBeamSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4CEmcTestBeamSteppingAction( PHG4CEmcTestBeamDetector* );

  //! destroctor
  virtual ~PHG4CEmcTestBeamSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  private:

  //! pointer to the detector
  PHG4CEmcTestBeamDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4HitContainer * absorberhits_;
  PHG4Hit *hit;
};


#endif // PHG4CEmcTestBeamSteppingAction_h
