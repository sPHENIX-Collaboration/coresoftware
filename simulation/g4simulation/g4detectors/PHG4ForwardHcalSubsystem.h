// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4FORWARDHCALSUBSYSTEM_H
#define G4DETECTORS_PHG4FORWARDHCALSUBSYSTEM_H

#include <g4main/PHG4Subsystem.h>

#include <string>                  // for string

class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4ForwardHcalDetector;
class PHG4SteppingAction;

class PHG4ForwardHcalSubsystem : public PHG4Subsystem
{
 public:
  /** Constructor
   */
  PHG4ForwardHcalSubsystem(const std::string &name = "FORWARD_HCAL_DEFAULT", const int layer = 0);

  /** Destructor
   */
  virtual ~PHG4ForwardHcalSubsystem();

  /**
     Creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
     Creates the stepping action and place it on the node tree, under "ACTIONS" node
     Creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int Init(PHCompositeNode *);

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
  void SetTowerMappingFile(const std::string &filename)
  {
    mappingfile_ = filename;
  }

  void SetActive(const int i = 1) { active = i; }
  void SetAbsorberActive(const int i = 1) { absorber_active = i; }
  void BlackHole(const int i = 1) { blackhole = i; }

 private:
  /** Pointer to the Geant4 implementation of the detector
   */
  PHG4ForwardHcalDetector *m_Detector;

  /** Stepping action
   */
  PHG4SteppingAction *m_SteppingAction;
  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction *m_DisplayAction;

  int active;
  int absorber_active;
  int blackhole;

  std::string detector_type;
  std::string mappingfile_;
};

#endif
