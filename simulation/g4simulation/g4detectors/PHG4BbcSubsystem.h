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

#ifndef G4DETECTORS_PHG4BBCSUBSYSTEM_H
#define G4DETECTORS_PHG4BBCSUBSYSTEM_H

#include "PHG4DetectorSubsystem.h"

#include <string>

class PHCompositeNode;
class PHG4Detector;
class PHG4BbcDetector;
class PHG4SteppingAction;

/**
   * \brief Fun4All module to simulate the BBC detector, aka MBD.
   *
   * The detector is constructed and registered via PHG4BbcDetector
   *
   * The PHG4SteppingAction needs to be updated more, but will provide the info for the hit time
   *
   * \see PHG4BbcDetector
   * \see PHG4BbcSubsystem
   *
   */
class PHG4BbcSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4BbcSubsystem(const std::string& name = "BBC");

  //! destructor
  ~PHG4BbcSubsystem(void) override
  {
  }

  /*!
  creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
  creates the stepping action and place it on the node tree, under "ACTIONS" node
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
  PHG4Detector* GetDetector() const override;
  PHG4SteppingAction* GetSteppingAction(void) const override;

  //! Print info (from SubsysReco)
  void Print(const std::string& what = "ALL") const override;

 protected:
  // \brief Set default parameter values
  void SetDefaultParameters() override;

  //! MBD geometry and construction
  /*! derives from PHG4Detector */
  PHG4BbcDetector* m_detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* m_steppingAction;

};

#endif
