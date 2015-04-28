#ifndef PHG4VSiliconTrackerSteppingAction_h
#define PHG4VSiliconTrackerSteppingAction_h

#include "g4main/PHG4SteppingAction.h"

class PHG4SiliconTrackerDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4SiliconTrackerSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4SiliconTrackerSteppingAction( PHG4SiliconTrackerDetector* );

  //! destroctor
  virtual ~PHG4SiliconTrackerSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  private:

  //! pointer to the detector
  PHG4SiliconTrackerDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4HitContainer * absorberhits_;
  PHG4Hit *hit;
};


#endif // PHG4SiliconTrackerSteppingAction_h
