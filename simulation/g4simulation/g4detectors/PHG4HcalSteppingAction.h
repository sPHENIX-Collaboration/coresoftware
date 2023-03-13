// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4HCALSTEPPINGACTION_H
#define G4DETECTORS_PHG4HCALSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <cmath>

class G4Step;
class PHCompositeNode;
class PHG4HcalDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;

class PHG4HcalSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  explicit PHG4HcalSteppingAction(PHG4HcalDetector*);

  //! destroctor
  ~PHG4HcalSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step*, bool) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode*) override;

  void set_zmin(const float z) { zmin = z; }
  void set_zmax(const float z) { zmax = z; }

  void SetLightScintModel(const bool b = true)
  {
    light_scint_model_ = b;
  }

 private:
  //! pointer to the detector
  PHG4HcalDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer* m_HitContainer = nullptr;
  PHG4HitContainer* m_AbsorberHits = nullptr;
  PHG4HitContainer* m_SaveHitContainer = nullptr;
  PHG4Hit* m_Hit = nullptr;
  PHG4Shower* m_SaveShower = nullptr;
  float zmin = NAN;
  float zmax = NAN;

  bool light_scint_model_ = true;
};

#endif  // PHG4VHcalSteppingAction_h
