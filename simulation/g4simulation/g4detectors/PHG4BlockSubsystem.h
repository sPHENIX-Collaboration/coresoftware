#ifndef PHG4BlockSubsystem_h
#define PHG4BlockSubsystem_h

#include "PHG4DetectorSubsystem.h"

class PHG4BlockDetector;
class PHG4BlockSteppingAction;

class PHG4BlockSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4BlockSubsystem(const std::string& name = "BLOCK", const int layer = 0);

  //! destructor
  virtual ~PHG4BlockSubsystem(void)
  {
  }

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
  PHG4Detector* GetDetector(void) const;
  PHG4SteppingAction* GetSteppingAction(void) const { return _steppingAction; }
 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4BlockDetector* _detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* _steppingAction;
};

#endif
