// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BEAMLINEMAGNETSUBSYSTEM_H
#define G4DETECTORS_PHG4BEAMLINEMAGNETSUBSYSTEM_H

#include "PHG4DetectorSubsystem.h"

#include <string>  // for string

class PHCompositeNode;
class PHG4BeamlineMagnetDetector;
class PHG4Detector;

class PHG4BeamlineMagnetSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4BeamlineMagnetSubsystem(const std::string &name = "CYLINDER", const int layer = 0);

  //! destructor
  ~PHG4BeamlineMagnetSubsystem(void) override
  {
  }

  //! init runwise stuff
  /*!
  creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode *) override;

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode *) override;

  //! Print info (from SubsysReco)
  void Print(const std::string &what = "ALL") const override;

  //! accessors (reimplemented)
  PHG4Detector *GetDetector(void) const override;

 private:
  void SetDefaultParameters() override;

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4BeamlineMagnetDetector *detector_;
};

#endif
