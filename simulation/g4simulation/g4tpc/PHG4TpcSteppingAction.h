#ifndef G4TPC_PHG4VTPCSTEPPINGACTION_H
#define G4TPC_PHG4VTPCSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4VPhysicalVolume;
class PHG4TpcDetector;
class PHParameters;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;

class PHG4TpcSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4TpcSteppingAction(PHG4TpcDetector *, const PHParameters *parameters);

  //! destructor
  virtual ~PHG4TpcSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step *, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode *);

 private:
  //! pointer to the detector
  PHG4TpcDetector *detector_;

  //! pointer to hit container
  PHG4HitContainer *hits_;
  PHG4HitContainer *absorberhits_;
  PHG4Hit *hit;
  const PHParameters *params;
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

#endif  // G4TPC_PHG4TPCSTEPPINGACTION_H
