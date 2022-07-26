// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4ENVELOPESTEPPINGACTION_H
#define G4DETECTORS_PHG4ENVELOPESTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class PHCompositeNode;
class PHG4EnvelopeDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4EnvelopeSteppingAction : public PHG4SteppingAction
{
 public:
  //Constructor
  explicit PHG4EnvelopeSteppingAction(PHG4EnvelopeDetector*);

  //Destructor
  ~PHG4EnvelopeSteppingAction() override
  {
  }

  //Stepping Action
  bool UserSteppingAction(const G4Step*, bool) override;

  //reimplemented from base class
  void SetInterfacePointers(PHCompositeNode*) override;

 private:
  //pointer to the detector
  PHG4EnvelopeDetector* detector_;

  //pointer to hit container
  PHG4HitContainer* hits_;
  PHG4Hit* hit;
};

#endif  //PHG4EnvelopeSteppingAction_h
