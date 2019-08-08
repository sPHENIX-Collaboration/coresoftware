// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MVTX_PHG4OUTERTRACKERSUBSYSTEM_H
#define G4MVTX_PHG4OUTERTRACKERSUBSYSTEM_H

#include <g4detectors/PHG4DetectorGroupSubsystem.h>

#include <string>                                    // for string

class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4OuterTrackerDetector;
class PHG4SteppingAction;

class PHG4OuterTrackerSubsystem : public PHG4DetectorGroupSubsystem
{
 public:
  //! constructor
  PHG4OuterTrackerSubsystem(const std::string& name = "PHG4OuterTrackerSubsystem", const int layer = 0);

  //! destructor
  virtual ~PHG4OuterTrackerSubsystem();

  //! InitRunSubsystem
  /*!
  called during InitRun (the original InitRun does common setup and calls this one)
  creates the detector object 
  ceates the stepping action 
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
  virtual PHG4Detector* GetDetector(void) const;
  virtual PHG4SteppingAction* GetSteppingAction(void) const;

  PHG4DisplayAction* GetDisplayAction() const { return m_DisplayAction; }

 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4OuterTrackerDetector* m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* steppingAction_;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction;

  // These are passed on to the detector class
  unsigned int layer;

  std::string detector_type;
};

#endif
