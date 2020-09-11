// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/* vim: set sw=2 ft=cpp: */

#ifndef G4EPD_PHG4EPSTEPPINGACTION_H
#define G4EPD_PHG4EPSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <cstdint>

class G4Step;
class PHCompositeNode;
class PHG4EPDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHParametersContainer;

class PHG4EPSteppingAction : public PHG4SteppingAction
{
 public:
  PHG4EPSteppingAction(PHG4EPDetector*, PHParametersContainer const*);
  ~PHG4EPSteppingAction() override;

  bool UserSteppingAction(const G4Step*, bool) override;

  void SetInterfacePointers(PHCompositeNode*) override;

 private:
  PHG4EPDetector* m_detector;

  PHG4HitContainer* m_hit_container;
  PHG4Hit* m_hit;

  int32_t poststatus;
};

#endif /* G4EPD_PHG4EPSTEPPINGACTION_H */
