// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/* vim: set sw=2 ft=cpp: */

#ifndef G4EPD_PHG4EPDSUBSYSTEM_H
#define G4EPD_PHG4EPDSUBSYSTEM_H

#include <g4detectors/PHG4DetectorGroupSubsystem.h>

#include <cstdint>
#include <string>

class PHCompositeNode;
class PHG4Detector;
class PHG4EPDetector;
class PHG4EPSteppingAction;
class PHG4SteppingAction;

class PHG4EPDSubsystem : public PHG4DetectorGroupSubsystem
{
 public:
  PHG4EPDSubsystem(std::string const& name);

  int32_t InitRunSubsystem(PHCompositeNode* node) override;

  int32_t process_event(PHCompositeNode* node) override;

  PHG4Detector* GetDetector() const override;
  PHG4SteppingAction* GetSteppingAction() const override;

 private:
  void SetDefaultParameters() override;

  PHG4EPDetector* m_detector;
  PHG4EPSteppingAction* m_stepaction;
};

#endif /* G4EPD_PHG4EPDSUBSYSTEM_H */
