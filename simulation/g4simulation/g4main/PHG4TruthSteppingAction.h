#ifndef PHG4TruthSteppingAction_h
#define PHG4TruthSteppingAction_h

#include "PHG4SteppingAction.h"

class PHG4TruthEventAction;

class PHG4TruthSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4TruthSteppingAction( PHG4TruthEventAction* );

  //! destroctor
  virtual ~PHG4TruthSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool );

  private:

  //! pointer to the detector
  PHG4TruthEventAction* eventAction_;

};


#endif
