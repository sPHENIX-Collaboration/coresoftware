// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/* vim: set sw=2 ft=cpp: */

#ifndef G4EPD_PHG4EPDSUBSYSTEM_H
#define G4EPD_PHG4EPDSUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

#include <string>

class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4EPDDetector;
class PHG4SteppingAction;

class PHG4EPDSubsystem : public PHG4DetectorSubsystem
{
 public:
  PHG4EPDSubsystem(const std::string& name = "EPD");
  ~PHG4EPDSubsystem() override;

  /**
     Creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode* node) override;

  /** Event processing
   */
  int process_event(PHCompositeNode*) override;

  PHG4Detector* GetDetector() const override;
  PHG4SteppingAction* GetSteppingAction() const override { return m_SteppingAction; };
  PHG4DisplayAction* GetDisplayAction() const override { return m_DisplayAction; }

 private:
  void SetDefaultParameters() override;

  /** Pointer to the Geant4 implementation of the detector
   */
  PHG4EPDDetector* m_Detector = nullptr;

  /** Stepping action
   */
  PHG4SteppingAction* m_SteppingAction = nullptr;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction = nullptr;

  std::string m_HitNodeName;
  std::string m_SupportNodeName;
};

#endif /* G4EPD_PHG4EPDSUBSYSTEM_H */
