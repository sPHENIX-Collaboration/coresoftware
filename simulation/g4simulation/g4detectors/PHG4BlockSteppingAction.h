#ifndef PHG4VBlockSteppingAction_h
#define PHG4VBlockSteppingAction_h

#include "g4main/PHG4SteppingAction.h"

class PHG4BlockDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Parameters;

class PHG4BlockSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4BlockSteppingAction( PHG4BlockDetector*, const PHG4Parameters *parameters );

  //! destroctor
  virtual ~PHG4BlockSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  private:

  //! pointer to the detector
  PHG4BlockDetector* detector_;
  const PHG4Parameters *params;
  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4Hit *hit;

  int absorbertruth;
  int active;
  int IsBlackHole;
  
  int use_g4_steps;
};


#endif //__G4PHPHYTHIAREADER_H__
