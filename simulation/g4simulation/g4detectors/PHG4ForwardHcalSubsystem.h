// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4FORWARDHCALSUBSYSTEM_H
#define G4DETECTORS_PHG4FORWARDHCALSUBSYSTEM_H

#include "PHG4DetectorSubsystem.h"

#include <string>                  // for string

class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4ForwardHcalDetector;
class PHG4SteppingAction;

class PHG4ForwardHcalSubsystem : public PHG4DetectorSubsystem
{
 public:
  /** Constructor
   */
  PHG4ForwardHcalSubsystem(const std::string &name = "FORWARD_HCAL_DEFAULT", const int layer = 0);

  /** Destructor
   */
  virtual ~PHG4ForwardHcalSubsystem();

  /**
     Creates the m_Detector object
     Creates the stepping action
     Creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode *);

  /** Event processing
   */
  int process_event(PHCompositeNode *);

  /** Accessors (reimplemented)
   */
  PHG4Detector *GetDetector() const;
  PHG4SteppingAction *GetSteppingAction() const { return m_SteppingAction; }
  PHG4DisplayAction *GetDisplayAction() const { return m_DisplayAction; }

  /** Set mapping file for calorimeter towers
   */
  void SetTowerMappingFile(const std::string &filename);

 private:
  void SetDefaultParameters();

  /** Pointer to the Geant4 implementation of the detector
   */
  PHG4ForwardHcalDetector *m_Detector = nullptr;

  /** Stepping action
   */
  PHG4SteppingAction *m_SteppingAction = nullptr;
  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction *m_DisplayAction = nullptr;

  int blackhole;

  std::string detector_type;
  std::string mappingfile_;
};

#endif
