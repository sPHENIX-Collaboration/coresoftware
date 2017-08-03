#ifndef __PHG4PSTOFSTEPPINGACTION_H__
#define __PHG4PSTOFSTEPPINGACTION_H__

#include "g4main/PHG4SteppingAction.h"

class PHG4PSTOFDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4ParametersContainer;

class PHG4PSTOFSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4PSTOFSteppingAction(PHG4PSTOFDetector*, const PHG4ParametersContainer*);

  //! destroctor
  virtual ~PHG4PSTOFSteppingAction()
  {
  }

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*);

 private:
  //! pointer to the detector
  PHG4PSTOFDetector* detector_;
  const PHG4ParametersContainer* params;
  //! pointer to hit container
  PHG4HitContainer* hits_;
  PHG4Hit* hit;

  int active;
  int use_g4_steps;
};

#endif  // __PHG4PSTOFSTEPPINGACTION_H__
