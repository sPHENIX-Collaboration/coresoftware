// $Id$

/*!
 * \file ${file_name}
 * \brief
 * \author Mickey Chiu <chiu@bnl.gov>
 * \version $Revision$
 * \date $Date$
 */

#ifndef PHG4PSTOFSubsystem_h
#define PHG4PSTOFSubsystem_h

#include "PHG4DetectorGroupSubsystem.h"

#include <string>

class PHG4PSTOFDetector;

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
  virtual ~PHG4PSTOFSubsystem(void)
  {
  }

  /*!
  creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  virtual int InitRunSubsystem(PHCompositeNode*);

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  virtual int process_event(PHCompositeNode*);

  //! accessors (reimplemented)
  virtual PHG4Detector* GetDetector(void) const;
  virtual PHG4SteppingAction* GetSteppingAction(void) const;
  //! Print info (from SubsysReco)
  virtual void Print(const std::string& what = "ALL") const;

 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4PSTOFDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* steppingAction_;

};

#endif
