// This is the steppingaction heaer file for the hcal prototype
// created on 1/27/2014, Liang, HeXC
//
#ifndef G4DETECTORS_PHG4HCALPROTOTYPESTEPPINGACTION_H
#define G4DETECTORS_PHG4HCALPROTOTYPESTEPPINGACTION_H

#include "g4main/PHG4SteppingAction.h"

class G4Step;
class PHCompositeNode;
class PHG4HcalPrototypeDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4HcalPrototypeSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4HcalPrototypeSteppingAction( PHG4HcalPrototypeDetector* );

  //! destroctor
  virtual ~PHG4HcalPrototypeSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  private:

  //! pointer to the detector
  PHG4HcalPrototypeDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4HitContainer * absorberhits_;
  PHG4Hit *hit;
};


#endif // PHG4HcalPrototypeSteppingAction_h
