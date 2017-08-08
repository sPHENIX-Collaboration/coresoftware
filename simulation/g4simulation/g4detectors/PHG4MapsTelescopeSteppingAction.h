#ifndef PHG4VMapsSteppingAction_h
#define PHG4VMapsSteppingAction_h

#include "g4main/PHG4SteppingAction.h"

class PHG4MapsTelescopeDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4MapsTelescopeSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4MapsTelescopeSteppingAction( PHG4MapsTelescopeDetector* );

  //! destroctor
  virtual ~PHG4MapsTelescopeSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  private:

  //! pointer to the detector
  PHG4MapsTelescopeDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4HitContainer * absorberhits_;
  PHG4Hit *hit;
};


#endif // PHG4MapsTelescopeSteppingAction_h
