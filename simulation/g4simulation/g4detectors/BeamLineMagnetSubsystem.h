// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_BEAMLINEMAGNETSUBSYSTEM_H
#define G4DETECTORS_BEAMLINEMAGNETSUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

#include <string>  // for string

class BeamLineMagnetDetector;
class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4SteppingAction;

class BeamLineMagnetSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  BeamLineMagnetSubsystem(const std::string& name = "CYLINDER", const int layer = 0);

  //! destructor
  virtual ~BeamLineMagnetSubsystem();

  //! init runwise stuff
  /*!
  creates the m_Detector object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode*);

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode*);

  //! Print info (from SubsysReco)
  void Print(const std::string& what = "ALL") const;

  //! accessors (reimplemented)
  PHG4Detector* GetDetector(void) const;
  PHG4SteppingAction* GetSteppingAction(void) const { return m_SteppingAction; }
  PHG4DisplayAction* GetDisplayAction() const { return m_DisplayAction; }

 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! defives from PHG4Detector */
  BeamLineMagnetDetector* m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* m_SteppingAction;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction;
};

#endif  // G4DETECTORS_BEAMLINEMAGNETSUBSYSTEM_H
