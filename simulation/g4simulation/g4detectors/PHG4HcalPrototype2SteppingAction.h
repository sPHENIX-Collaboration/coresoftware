// This is the steppingaction heaer file for the hcal prototype
// created on 1/27/2014, Liang, HeXC
//
#ifndef PHG4VHcalPrototype2SteppingAction_h
#define PHG4VHcalPrototype2SteppingAction_h

#include "g4main/PHG4SteppingAction.h"

class PHG4HcalPrototype2Detector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4HcalPrototype2SteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4HcalPrototype2SteppingAction( PHG4HcalPrototype2Detector* );

  //! destroctor
  virtual ~PHG4HcalPrototype2SteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  private:

  //! pointer to the detector
  PHG4HcalPrototype2Detector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4HitContainer * absorberhits_;
  PHG4Hit *hit;
};


#endif // PHG4HcalPrototype2SteppingAction_h
