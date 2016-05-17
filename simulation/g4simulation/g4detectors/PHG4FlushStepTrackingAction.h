#ifndef PHG4FlushStepTrackingAction_h__
#define PHG4FlushStepTrackingAction_h__

#include <g4main/PHG4TrackingAction.h>

class G4Track;
class PHG4SteppingAction;

class PHG4FlushStepTrackingAction: public PHG4TrackingAction
{
public:
  PHG4FlushStepTrackingAction (PHG4SteppingAction *);

  virtual ~PHG4FlushStepTrackingAction() {}

  void PreUserTrackingAction(const G4Track*) {return;}

  void PostUserTrackingAction(const G4Track*);

 private:
  PHG4SteppingAction *steppingAction_;
};

#endif
