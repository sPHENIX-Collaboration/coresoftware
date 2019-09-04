// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4HCALSTEPPINGACTION_H
#define G4DETECTORS_PHG4HCALSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class PHCompositeNode;
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
  PHG4HitContainer *m_SaveHitContainer;
  PHG4Hit *hit;
  float zmin;
  float zmax;

  bool  light_scint_model_;
};


#endif // PHG4VHcalSteppingAction_h
