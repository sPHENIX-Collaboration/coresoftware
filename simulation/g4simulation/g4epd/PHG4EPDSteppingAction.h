// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/* vim: set sw=2 ft=cpp: */

#ifndef G4EPD_PHG4EPDSTEPPINGACTION_H
#define G4EPD_PHG4EPDSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <cstdint>

class G4Step;
class PHCompositeNode;
class PHG4EPDDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHParametersContainer;

class PHG4EPDSteppingAction : public PHG4SteppingAction
{
 public:
  PHG4EPDSteppingAction(PHG4EPDDetector*, PHParametersContainer const*);
  ~PHG4EPDSteppingAction() override;

  bool UserSteppingAction(const G4Step*, bool) override;

  void SetInterfacePointers(PHCompositeNode*) override;

 private:
  PHG4EPDDetector* m_detector;

  PHG4HitContainer* m_hit_container;
  PHG4Hit* m_hit;

  int32_t poststatus;
};

#endif /* G4EPD_PHG4EPSTEPPINGACTION_H */
