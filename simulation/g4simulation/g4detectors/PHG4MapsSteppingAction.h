#ifndef PHG4VMapsSteppingAction_h
#define PHG4VMapsSteppingAction_h

#include <g4main/PHG4SteppingAction.h>

class PHG4MapsDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4MapsSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4MapsSteppingAction( PHG4MapsDetector* );

  //! destroctor
  virtual ~PHG4MapsSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  private:

  //! pointer to the detector
  PHG4MapsDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4HitContainer * absorberhits_;
  PHG4Hit *hit;
};


#endif // PHG4MapsSteppingAction_h
