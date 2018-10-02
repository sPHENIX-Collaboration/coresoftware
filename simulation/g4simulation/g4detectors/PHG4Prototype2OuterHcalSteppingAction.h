#ifndef PHG4VPrototype2OuterHcalSteppingAction_h
#define PHG4VPrototype2OuterHcalSteppingAction_h

#include <g4main/PHG4SteppingAction.h>

class PHG4Prototype2OuterHcalDetector;
class PHParameters;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;

class PHG4Prototype2OuterHcalSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4Prototype2OuterHcalSteppingAction(PHG4Prototype2OuterHcalDetector *, const PHParameters *parameters);

  //! dtor
  virtual ~PHG4Prototype2OuterHcalSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step *, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode *);

  double GetLightCorrection(const double r) const;

 private:
  //! pointer to the detector
  PHG4Prototype2OuterHcalDetector *m_Detector;

  //! pointer to hit container
  PHG4HitContainer *hits_;
  PHG4HitContainer *absorberhits_;
  PHG4Hit *hit;
  const PHParameters *params;
  PHG4HitContainer *savehitcontainer;
  PHG4Shower *saveshower;
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

#endif  // PHG4Prototype2OuterHcalSteppingAction_h
