#ifndef PHG4CylinderSubsystem_h
#define PHG4CylinderSubsystem_h

#include "PHG4DetectorSubsystem.h"

#include <Geant4/G4String.hh>
#include <Geant4/G4Types.hh>

class PHG4CylinderDetector;
class PHG4SteppingAction;

class PHG4CylinderSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4CylinderSubsystem(const std::string& name = "CYLINDER", const int layer = 0);

  //! destructor
  virtual ~PHG4CylinderSubsystem(void)
  {
  }

  //! init runwise stuff
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

  //! Print info (from SubsysReco)
  void Print(const std::string& what = "ALL") const;

  //! accessors (reimplemented)
  PHG4Detector* GetDetector(void) const;
  PHG4SteppingAction* GetSteppingAction(void) const { return steppingAction_; }
 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4CylinderDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* steppingAction_;
};

#endif
