#ifndef __TPCSTEPPINGACTION_H__
#define __TPCSTEPPINGACTION_H__

#include <g4main/PHG4SteppingAction.h>

class TPCDetector;
class TPCHitsContainer;
class TPCHit;

class TPCSteppingAction : public PHG4SteppingAction {
 public:
  TPCSteppingAction( TPCDetector* );
  virtual ~TPCSteppingAction();
  virtual bool UserSteppingAction(const G4Step*, bool);
  virtual void SetInterfacePointers(PHCompositeNode*);
  void SetVerbosity(int v) {fVerbosity=v;}

 private:
  int fVerbosity;
  TPCDetector *fDetector;
  TPCHitsContainer *fHits;
  TPCHit *fHit;
};

#endif
