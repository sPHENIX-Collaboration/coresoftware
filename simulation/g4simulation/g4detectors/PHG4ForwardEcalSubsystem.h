// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4FORWARDECALSUBSYSTEM_H
#define G4DETECTORS_PHG4FORWARDECALSUBSYSTEM_H

#include "PHG4DetectorSubsystem.h"

#include <string>  // for string

class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4ForwardEcalDetector;
class PHG4SteppingAction;

class PHG4ForwardEcalSubsystem : public PHG4DetectorSubsystem
{
 public:
  /** Constructor
   */
  PHG4ForwardEcalSubsystem(const std::string& name = "FORWARD_ECAL_DEFAULT", const int layer = 0);

  /** Destructor
   */
  virtual ~PHG4ForwardEcalSubsystem();

  /**
     Creates the m_Detector object and place it on the node tree, under "DETECTORS" node (or whatever)
     Creates the stepping action and place it on the node tree, under "ACTIONS" node
     Creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode*);

  /** Event processing
   */
  int process_event(PHCompositeNode*);

  /** Accessors (reimplemented)
   */
  PHG4Detector* GetDetector() const;
  PHG4SteppingAction* GetSteppingAction() const { return m_SteppingAction; }

  PHG4DisplayAction* GetDisplayAction() const { return m_DisplayAction; }

  void SetEICDetector() { m_EICDetectorFlag = 1; }
  void SetfsPHENIXDetector() { m_EICDetectorFlag = 0; }
  void SetTowerMappingFile(const std::string& filename);

 private:
  void SetDefaultParameters();

  /** Pointer to the Geant4 implementation of the detector
   */
  PHG4ForwardEcalDetector* m_Detector;

  /** Stepping action
   */
  PHG4SteppingAction* m_SteppingAction;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction;

  int m_EICDetectorFlag;
};

#endif
