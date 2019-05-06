// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4MVTX_PHG4VMVTXSTEPPINGACTION_H
#define G4MVTX_PHG4VMVTXSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class PHG4MVTXDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4MVTXSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4MVTXSteppingAction(PHG4MVTXDetector*);

  //! destroctor
  virtual ~PHG4MVTXSteppingAction()
  {
  }

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*);

 private:
  //! pointer to the detector
  PHG4MVTXDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer* hits_;
  PHG4HitContainer* absorberhits_;
  PHG4Hit* hit;
};

#endif  // G4MVTX_PHG4VMVTXSTEPPINGACTION_H
