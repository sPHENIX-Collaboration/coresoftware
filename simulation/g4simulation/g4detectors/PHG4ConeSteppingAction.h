// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CONESTEPPINGACTION_H
#define G4DETECTORS_PHG4CONESTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class PHCompositeNode;
class PHG4ConeDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;

class PHG4ConeSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  explicit PHG4ConeSteppingAction(PHG4ConeDetector*);

  //! destructor
  ~PHG4ConeSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step*, bool) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode*) override;

 private:
  //! pointer to the detector
  PHG4ConeDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer* hits_;
  PHG4Hit* hit;
  PHG4Shower* saveshower;
};

#endif  // G4DETECTORS_PHG4CONESTEPPINGACTION_H
