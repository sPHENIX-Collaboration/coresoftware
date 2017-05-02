#ifndef __TPCSTEPPINGACTION_H__
#define __TPCSTEPPINGACTION_H__

#include <g4main/PHG4SteppingAction.h>

class TPCDetector;
class TPCHitsContainer;
class TPCHit;
class PHG4HitContainer;

class TPCSteppingAction : public PHG4SteppingAction {
 public:
  TPCSteppingAction( TPCDetector* );
  virtual ~TPCSteppingAction();
  virtual bool UserSteppingAction(const G4Step*, bool);
  virtual void SetInterfacePointers(PHCompositeNode*);
  void SetVerbosity(int v) {fVerbosity=v;}

 private:
  bool SteppingActionTPC(const G4Step*);
  void SteppingActionPHG4(const G4Step*);

  int fVerbosity;
  TPCDetector *fDetector;
  TPCHitsContainer *fHits;
  TPCHit *fHit;
  PHG4HitContainer *fPHHits;
};

#endif
