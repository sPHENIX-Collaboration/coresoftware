#ifndef PHG4VOuterHcalPrototype2SteppingAction_h
#define PHG4VOuterHcalPrototype2SteppingAction_h

#include <g4main/PHG4SteppingAction.h>

class PHG4OuterHcalPrototype2Detector;
class PHG4Parameters;
class PHG4Hit;
class PHG4HitContainer;

class PHG4OuterHcalPrototype2SteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4OuterHcalPrototype2SteppingAction( PHG4OuterHcalPrototype2Detector*, PHG4Parameters *parameters );

  //! destroctor
  virtual ~PHG4OuterHcalPrototype2SteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  double GetLightCorrection(const double r) const;

  private:

  //! pointer to the detector
  PHG4OuterHcalPrototype2Detector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4HitContainer * absorberhits_;
  PHG4Hit *hit;
  PHG4Parameters *params;
  // since getting parameters is a map search we do not want to
  // do this in every step, the parameters used are cached
  // in the following variables
  int absorbertruth;
  int IsActive;
  int IsBlackHole;
  int light_scint_model;
  
  double light_balance_inner_corr;
  double light_balance_inner_radius;
  double light_balance_outer_corr;
  double light_balance_outer_radius;
};


#endif // PHG4OuterHcalPrototype2SteppingAction_h
