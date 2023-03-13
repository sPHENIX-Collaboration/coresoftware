// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CEMCTESTBEAMSTEPPINGACTION_H
#define G4DETECTORS_PHG4CEMCTESTBEAMSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class PHCompositeNode;
class PHG4CEmcTestBeamDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4CEmcTestBeamSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  explicit PHG4CEmcTestBeamSteppingAction(PHG4CEmcTestBeamDetector*);

  //! destroctor
  ~PHG4CEmcTestBeamSteppingAction() override
  {
  }

  //! stepping action
  bool UserSteppingAction(const G4Step*, bool) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode*) override;

 private:
  //! pointer to the detector
  PHG4CEmcTestBeamDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer* hits_;
  PHG4HitContainer* absorberhits_;
  PHG4Hit* hit;
};

#endif  // PHG4CEmcTestBeamSteppingAction_h
