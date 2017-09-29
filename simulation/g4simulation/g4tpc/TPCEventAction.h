#ifndef __TPCEVENTACTION_H__
#define __TPCEVENTACTION_H__

#include <g4main/PHG4EventAction.h>

class G4Event;
class PHCompositeNode;

class TPCEventAction : public PHG4EventAction {
 public:
  TPCEventAction(PHCompositeNode *topNode);
  virtual ~TPCEventAction() {}
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
  void SetVerbosity(int v) {fVerbosity=v;}
  
 private:
  int fVerbosity;
  PHCompositeNode *fNode;
  
};


#endif
