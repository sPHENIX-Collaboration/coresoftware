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
class PHParametersContainer;

class PHG4MvtxSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4MvtxSteppingAction(PHG4MvtxDetector *detector, PHParametersContainer *paramscont);

  //! destroctor
  ~PHG4MvtxSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step *, bool) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode *) override;

  void SetHitNodeName(const std::string &type, const std::string &name) override;

 private:
  //! pointer to the detector
  PHG4MvtxDetector *m_Detector{nullptr};

  //! pointer to hit container
  PHG4HitContainer *m_HitContainer{nullptr};
  PHG4HitContainer *m_SupportHitContainer{nullptr};
  PHG4HitContainer *m_SaveHitContainer{nullptr};
  PHG4Hit *m_Hit{nullptr};
  PHG4Shower *m_SaveShower{nullptr};

  int m_SupportActiveFlag{0};
  int m_BlackHoleFlag{0};

  std::string m_HitNodeName;
  std::string m_SupportNodeName;
};

#endif  // G4MVTX_PHG4VMVTXSTEPPINGACTION_H
