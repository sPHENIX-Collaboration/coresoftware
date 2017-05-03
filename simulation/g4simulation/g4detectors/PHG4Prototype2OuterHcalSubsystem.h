#ifndef PHG4Prototype2OuterHcalSubsystem_h
#define PHG4Prototype2OuterHcalSubsystem_h

#include "PHG4DetectorSubsystem.h"

#include <map>
#include <set>
#include <string>

class PHG4Prototype2OuterHcalDetector;
class PHG4SteppingAction;

class PHG4Prototype2OuterHcalSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4Prototype2OuterHcalSubsystem(const std::string& name = "HCALIN", const int layer = 0);

  //! destructor
  virtual ~PHG4Prototype2OuterHcalSubsystem(void)
  {
  }

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
  virtual PHG4Detector* GetDetector(void) const;
  virtual PHG4SteppingAction* GetSteppingAction(void) const { return steppingAction_; }
  void SetLightCorrection(const double inner_radius, const double inner_corr, const double outer_radius, const double outer_corr);

 protected:
  void SetDefaultParameters();

  //! detector geometry
  /*! derives from PHG4Detector */
  PHG4Prototype2OuterHcalDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingAction */
  PHG4SteppingAction* steppingAction_;
};

#endif
