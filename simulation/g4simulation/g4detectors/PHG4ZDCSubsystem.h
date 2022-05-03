// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4FORWARDECALSUBSYSTEM_H
#define G4DETECTORS_PHG4FORWARDECALSUBSYSTEM_H

#include "PHG4DetectorSubsystem.h"

#include <string>  // for string

class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4ZDCDetector;
class PHG4SteppingAction;

class PHG4ZDCSubsystem : public PHG4DetectorSubsystem
{
 public:
  /** Constructor
   */
  explicit PHG4ZDCSubsystem(const std::string& name, const int layer);

  /** Destructor
   */
  ~PHG4ZDCSubsystem() override;

  /**
     Creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode*) override;

  /** Event processing
   */
  int process_event(PHCompositeNode*) override;

  /** Accessors (reimplemented)
   */
  PHG4Detector* GetDetector() const override;
  PHG4SteppingAction* GetSteppingAction() const override { return m_SteppingAction; }

  PHG4DisplayAction* GetDisplayAction() const override { return m_DisplayAction; }

 private:
  void SetDefaultParameters() override;

  /** Pointer to the Geant4 implementation of the detector
   */
  PHG4ZDCDetector* m_Detector = nullptr;

  /** Stepping action
   */
  PHG4SteppingAction* m_SteppingAction = nullptr;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction = nullptr;

  std::string m_HitNodeName;
  std::string m_AbsorberNodeName;
  std::string m_SupportNodeName;
};

#endif
