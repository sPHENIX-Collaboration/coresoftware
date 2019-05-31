#ifndef G4DETECTORS_PHG4FCALSTEPPINGACTION_H
#define G4DETECTORS_PHG4FCALSTEPPINGACTION_H

#include <Geant4/G4UserSteppingAction.hh>

class PHCompositeNode;
class PHG4FCalDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4FCalSteppingAction : public G4UserSteppingAction
{
  public:
    PHG4FCalSteppingAction( PHG4FCalDetector* );
    virtual ~PHG4FCalSteppingAction(){}
    
    virtual void UserSteppingAction(const G4Step*);
    
    virtual void SetInterfacePointers( PHCompositeNode* );
    
  private:
    PHG4FCalDetector* detector_;
    PHG4HitContainer* hits_;
    PHG4Hit* hit;
};



#endif

