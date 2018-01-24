#ifndef PHG4VTPCSteppingAction_h
#define PHG4VTPCSteppingAction_h

#include <g4main/PHG4SteppingAction.h>

class G4VPhysicalVolume;
class PHG4TPCDetector;
class PHG4Parameters;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;

class PHG4TPCSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4TPCSteppingAction(PHG4TPCDetector *, const PHG4Parameters *parameters);

  //! destructor
  virtual ~PHG4TPCSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step *, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode *);

 private:
  //! pointer to the detector
  PHG4TPCDetector *detector_;

  //! pointer to hit container
  PHG4HitContainer *hits_;
  PHG4HitContainer *absorberhits_;
  PHG4Hit *hit;
  const PHG4Parameters *params;
  PHG4HitContainer *savehitcontainer;
  PHG4Shower *saveshower;
  G4VPhysicalVolume *savevolpre;
  G4VPhysicalVolume *savevolpost;
  int savetrackid;
  int saveprestepstatus;
  int savepoststepstatus;
  // since getting parameters is a map search we do not want to
  // do this in every step, the parameters used are cached
  // in the following variables
  int IsActive;
  int IsBlackHole;
  int use_g4_steps;
};

#endif  // PHG4TPCSteppingAction_h
