#ifndef PHG4VSiliconTrackerSteppingAction_h
#define PHG4VSiliconTrackerSteppingAction_h

#include <g4main/PHG4SteppingAction.h>

class PHG4SiliconTrackerDetector;
class PHG4Parameters;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;

class PHG4SiliconTrackerSteppingAction : public PHG4SteppingAction
{
 public:
  PHG4SiliconTrackerSteppingAction(PHG4SiliconTrackerDetector *, const PHG4Parameters *parameters);

  virtual ~PHG4SiliconTrackerSteppingAction();

  virtual bool UserSteppingAction(const G4Step *, bool);

  virtual void SetInterfacePointers(PHCompositeNode *);

 private:
  //! pointer to the detector
  PHG4SiliconTrackerDetector *detector_;

  //! pointer to hit container
  PHG4HitContainer *hits_;
  PHG4HitContainer *absorberhits_;
  PHG4Hit *hit;
  PHG4HitContainer *savehitcontainer;
  PHG4Shower *saveshower;
  const PHG4Parameters *params;

  int IsActive;
  int IsBlackHole;
  int toggle;
};

#endif  // PHG4SiliconTrackerSteppingAction_h
