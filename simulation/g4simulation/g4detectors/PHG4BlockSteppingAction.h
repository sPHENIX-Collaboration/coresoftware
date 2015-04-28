#ifndef PHG4VBlockSteppingAction_h
#define PHG4VBlockSteppingAction_h

#include "g4main/PHG4SteppingAction.h"

class PHG4BlockDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4BlockSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4BlockSteppingAction( PHG4BlockDetector* );

  //! destroctor
  virtual ~PHG4BlockSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );
  void UseG4Steps(const int i = 1) {use_g4_steps = i;}
  void UseIonizationEnergy(const int i) {use_ionisation_energy = i;}

  private:

  //! pointer to the detector
  PHG4BlockDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4Hit *hit;

  int use_g4_steps;
  int use_ionisation_energy;
};


#endif //__G4PHPHYTHIAREADER_H__
