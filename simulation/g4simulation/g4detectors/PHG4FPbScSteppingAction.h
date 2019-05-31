#ifndef G4DETECTORS_PHG4FPBSCSTEPPINGACTION_H
#define G4DETECTORS_PHG4FPBSCSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class PHCompositeNode;
class PHG4FPbScDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4FPbScSteppingAction : public PHG4SteppingAction
{
  public:
    PHG4FPbScSteppingAction( PHG4FPbScDetector* );
    virtual ~PHG4FPbScSteppingAction(){}
    
    virtual bool UserSteppingAction(const G4Step*, bool);
    
    virtual void SetInterfacePointers( PHCompositeNode* );
    
  private:
    PHG4FPbScDetector* detector_;
    PHG4HitContainer* hits_;
};



#endif

