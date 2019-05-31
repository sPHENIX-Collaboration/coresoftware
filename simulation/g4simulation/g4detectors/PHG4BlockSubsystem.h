// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BLOCKSUBSYSTEM_H
#define G4DETECTORS_PHG4BLOCKSUBSYSTEM_H

#include "PHG4DetectorSubsystem.h"

#include <string>                   // for string

class PHCompositeNode;
class PHG4Detector;
class PHG4BlockDetector;
class PHG4DisplayAction;
class PHG4SteppingAction;

class PHG4BlockSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4BlockSubsystem(const std::string& name = "BLOCK", const int layer = 0);

  //! destructor
  virtual ~PHG4BlockSubsystem();

  //! InitRunSubsystem
  /*!
  creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
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

  //! accessors (reimplemented)
  PHG4Detector* GetDetector() const;

  PHG4SteppingAction* GetSteppingAction() const { return m_SteppingAction; }

  PHG4DisplayAction* GetDisplayAction() const { return m_DisplayAction; }

 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! defines from PHG4Detector */
  PHG4BlockDetector* m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* m_SteppingAction;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction;
};

#endif
