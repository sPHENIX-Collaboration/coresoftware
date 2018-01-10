#ifndef PHG4TPCSubsystem_h
#define PHG4TPCSubsystem_h

#include <g4detectors/PHG4DetectorSubsystem.h>

#include <map>
#include <set>
#include <string>

class PHG4TPCDetector;
class PHG4Parameters;
class PHG4TPCSteppingAction;

class PHG4TPCSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4TPCSubsystem(const std::string &name = "TPC", const int layer = 0);

  //! destructor
  virtual ~PHG4TPCSubsystem(void)
  {
  }

  /*!
  creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
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
 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! derives from PHG4Detector */
  PHG4TPCDetector *detector_;

  //! detector "stepping" action, executes after every G4 step
  /*! derives from PHG4SteppingAction */
  PHG4SteppingAction *steppingAction_;

  //! detector event action executes before/after every event
  /*! derives from PHG4EventAction */
  PHG4EventAction *eventAction_;
};

#endif
