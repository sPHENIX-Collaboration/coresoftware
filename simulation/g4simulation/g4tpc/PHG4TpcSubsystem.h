// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4TPC_PHG4TPCSUBSYSTEM_H
#define G4TPC_PHG4TPCSUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

#include <string>

class PHG4DisplayAction;
class PHG4TpcDetector;
class PHG4TpcSteppingAction;

class PHG4TpcSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4TpcSubsystem(const std::string &name = "TPC", const int layer = 0);

  //! destructor
  virtual ~PHG4TpcSubsystem();

  /*!
  called during InitRun (the original InitRun does common setup and calls this one)
  creates the detector object 
  ceates the stepping action 
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode *);

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode *);

  //! Print info (from SubsysReco)
  void Print(const std::string &what = "ALL") const;

  //! accessors (reimplemented)
  PHG4Detector *GetDetector(void) const;

  PHG4SteppingAction *GetSteppingAction(void) const { return steppingAction_; }

  PHG4DisplayAction *GetDisplayAction() const { return m_DisplayAction; }

 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! derives from PHG4Detector */
  PHG4TpcDetector *detector_;

  //! detector "stepping" action, executes after every G4 step
  /*! derives from PHG4SteppingAction */
  PHG4SteppingAction *steppingAction_;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction *m_DisplayAction;
};

#endif
