#ifndef PHG4VInnerHcalSteppingAction_h
#define PHG4VInnerHcalSteppingAction_h

#include "g4main/PHG4SteppingAction.h"

class PHG4InnerHcalDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4InnerHcalSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4InnerHcalSteppingAction( PHG4InnerHcalDetector* );

  //! destroctor
  virtual ~PHG4InnerHcalSteppingAction()
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
  
  private:

  //! pointer to the detector
  PHG4InnerHcalDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4HitContainer * absorberhits_;
  PHG4Hit *hit;

  bool  light_balance_;
  float light_balance_inner_radius_;
  float light_balance_inner_corr_;
  float light_balance_outer_radius_;
  float light_balance_outer_corr_;
};


#endif // PHG4InnerHcalSteppingAction_h
