// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4MVTX_PHG4VMVTXSTEPPINGACTION_H
#define G4MVTX_PHG4VMVTXSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class PHCompositeNode;
class PHG4MvtxDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;

class PHG4MvtxSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4MvtxSteppingAction(PHG4MvtxDetector *);

  //! destroctor
  ~PHG4MvtxSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step *, bool) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode *) override;

 private:
  //! pointer to the detector
  PHG4MvtxDetector *m_Detector;

  //! pointer to hit container
  PHG4HitContainer *m_HitContainer;
  PHG4HitContainer *m_AbsorberhitContainer;
  PHG4Hit *m_Hit;
  PHG4Shower *m_SaveShower;
};

#endif  // G4MVTX_PHG4VMVTXSTEPPINGACTION_H
