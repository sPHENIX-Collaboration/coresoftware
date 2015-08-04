#ifndef PHG4VOuterHcalSteppingAction_h
#define PHG4VOuterHcalSteppingAction_h

#include "g4main/PHG4SteppingAction.h"

class PHG4OuterHcalDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4OuterHcalSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4OuterHcalSteppingAction( PHG4OuterHcalDetector* );

  //! destroctor
  virtual ~PHG4OuterHcalSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  float GetLightCorrection(float r);
  void SetLightCorrection(float inner_radius, float inner_corr,
			  float outer_radius, float outer_corr) {
    light_balance_ = true;
    light_balance_inner_radius_ = inner_radius;
    light_balance_inner_corr_ = inner_corr;
    light_balance_outer_radius_ = outer_radius;
    light_balance_outer_corr_ = outer_corr;
  }

  void SetLightScintModel(const bool b = true)
  {
    light_scint_model_ = b;
  }
  
  void FieldChecker (const G4Step*);

  private:

  //! pointer to the detector
  PHG4OuterHcalDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4HitContainer * absorberhits_;
  PHG4Hit *hit;

  bool  enable_field_checker_;
  bool  light_scint_model_;
  bool  light_balance_;
  float light_balance_inner_radius_;
  float light_balance_inner_corr_;
  float light_balance_outer_radius_;
  float light_balance_outer_corr_;
};


#endif // PHG4OuterHcalSteppingAction_h
