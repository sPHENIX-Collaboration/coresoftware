// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/* vim: set sw=2 ft=cpp: */

#ifndef G4EPD_PHG4EPDSTEPPINGACTION_H
#define G4EPD_PHG4EPDSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <Geant4/G4StepStatus.hh>

#include <string>

class G4Step;
class PHCompositeNode;
class PHG4EPDDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHParameters;

class PHG4EPDSteppingAction : public PHG4SteppingAction
{
 public:
  PHG4EPDSteppingAction(PHG4EPDDetector*, const PHParameters* parameters);
  ~PHG4EPDSteppingAction() override;

  bool UserSteppingAction(const G4Step*, bool) override;

  void SetInterfacePointers(PHCompositeNode*) override;

  void SetHitNodeName(const std::string& type, const std::string& name) override;

 private:
  PHG4EPDDetector* m_Detector = nullptr;

  PHG4HitContainer* m_HitContainer = nullptr;
  PHG4HitContainer* m_SupportHitContainer = nullptr;
  PHG4Hit* m_Hit = nullptr;

  G4StepStatus m_SavePostStepStatus = G4StepStatus::fUndefined;

  std::string m_HitNodeName;
  std::string m_SupportNodeName;
};

#endif /* G4EPD_PHG4EPSTEPPINGACTION_H */
