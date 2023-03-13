// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $Id$

/*!
 * \file ${file_name}
 * \brief
 * \author Mickey Chiu <chiu@bnl.gov>
 * \version $Revision$
 * \date $Date$
 */

#ifndef G4DETECTORS_PHG4PSTOFSUBSYSTEM_H
#define G4DETECTORS_PHG4PSTOFSUBSYSTEM_H

#include "PHG4DetectorGroupSubsystem.h"

#include <string>

class PHCompositeNode;
class PHG4Detector;
class PHG4PSTOFDetector;
class PHG4SteppingAction;

/**
   * \brief Fun4All module to simulate the Barrel PSTOF detector.
   *
   * The detector is constructed and registered via PHG4PSTOFDetector
   *
   * The PHG4SteppingAction needs to be written, but will provide the info for the hit time
   *
   * \see PHG4PSTOFDetector
   * \see PHG4PSTOFSubsystem
   *
   */
class PHG4PSTOFSubsystem : public PHG4DetectorGroupSubsystem
{
 public:
  //! constructor
  PHG4PSTOFSubsystem(const std::string& name = "PSTOF");

  //! destructor
  ~PHG4PSTOFSubsystem(void) override
  {
  }

  /*!
  creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
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

  //! accessors (reimplemented)
  PHG4Detector* GetDetector(void) const override;
  PHG4SteppingAction* GetSteppingAction(void) const override;
  //! Print info (from SubsysReco)
  void Print(const std::string& what = "ALL") const override;

 private:
  void SetDefaultParameters() override;

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4PSTOFDetector* detector_ = nullptr;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* steppingAction_ = nullptr;
};

#endif
