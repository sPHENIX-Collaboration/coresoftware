#ifndef G4DETECTORS_PHG4HCALSTEPPINGACTION_H
#define G4DETECTORS_PHG4HCALSTEPPINGACTION_H

#include "g4main/PHG4SteppingAction.h"
#include <string>

class PHG4HcalDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4HcalSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4HcalSteppingAction( PHG4HcalDetector* );

  //! destroctor
  virtual ~PHG4HcalSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  void set_zmin(const float z) {zmin = z;}
  void set_zmax(const float z) {zmax = z;}

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
  
  private:

  //! pointer to the detector
  PHG4HcalDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4HitContainer * absorberhits_;
  PHG4Hit *hit;
  float zmin;
  float zmax;

  bool  light_scint_model_;
  bool  light_balance_;
  float light_balance_inner_radius_;
  float light_balance_inner_corr_;
  float light_balance_outer_radius_;
  float light_balance_outer_corr_;
};


#endif // PHG4VHcalSteppingAction_h
