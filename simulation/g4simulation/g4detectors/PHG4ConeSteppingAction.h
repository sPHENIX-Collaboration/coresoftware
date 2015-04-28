#ifndef PHG4VConeSteppingAction_h
#define PHG4VConeSteppingAction_h

#include "g4main/PHG4SteppingAction.h"

class PHG4ConeDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4ConeSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4ConeSteppingAction( PHG4ConeDetector* );

  //! destroctor
  virtual ~PHG4ConeSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  private:

  //! pointer to the detector
  PHG4ConeDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4Hit *hit;
};


#endif //__G4PHPHYTHIAREADER_H__
