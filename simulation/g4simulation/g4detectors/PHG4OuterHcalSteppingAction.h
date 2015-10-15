#ifndef PHG4VOuterHcalSteppingAction_h
#define PHG4VOuterHcalSteppingAction_h

#include "g4main/PHG4SteppingAction.h"

class PHG4OuterHcalDetector;
class PHG4OuterHcalParameters;
class PHG4Hit;
class PHG4HitContainer;

class PHG4OuterHcalSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4OuterHcalSteppingAction( PHG4OuterHcalDetector* , PHG4OuterHcalParameters *parameters);

  //! destructor
  virtual ~PHG4OuterHcalSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  float GetLightCorrection(const float r) const;

  void FieldChecker (const G4Step*);

  private:

  //! pointer to the detector
  PHG4OuterHcalDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4HitContainer * absorberhits_;
  PHG4Hit *hit;

  PHG4OuterHcalParameters *params;
};


#endif // PHG4OuterHcalSteppingAction_h
