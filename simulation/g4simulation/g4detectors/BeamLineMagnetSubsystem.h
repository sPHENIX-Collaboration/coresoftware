// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_BEAMLINEMAGNETSUBSYSTEM_H
#define G4DETECTORS_BEAMLINEMAGNETSUBSYSTEM_H

#include "PHG4DetectorSubsystem.h"

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
  ~BeamLineMagnetSubsystem() override;

  //! init runwise stuff
  /*!
  creates the m_Detector object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode*) override;

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode*) override;

  //! Print info (from SubsysReco)
  void Print(const std::string& what = "ALL") const override;

  //! accessors (reimplemented)
  PHG4Detector* GetDetector(void) const override;
  PHG4SteppingAction* GetSteppingAction(void) const override { return m_SteppingAction; }
  PHG4DisplayAction* GetDisplayAction() const override { return m_DisplayAction; }

  // this method is used to check if it can be used as mothervolume
  // Subsystems which can be mothervolume need to implement this
  // and return true
  bool CanBeMotherSubsystem() const override { return true; }

 private:
  void SetDefaultParameters() override;

  //! detector geometry
  /*! defives from PHG4Detector */
  BeamLineMagnetDetector* m_Detector = nullptr;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* m_SteppingAction = nullptr;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction = nullptr;

  std::string m_HitNodeName;
  std::string m_AbsorberNodeName;
};

#endif  // G4DETECTORS_BEAMLINEMAGNETSUBSYSTEM_H
