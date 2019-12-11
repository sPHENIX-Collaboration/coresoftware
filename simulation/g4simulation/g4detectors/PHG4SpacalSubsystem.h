// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.2 $$
 * \date $$Date: 2014/08/12 03:49:12 $$
 */

#ifndef G4DETECTORS_PHG4SPACALSUBSYSTEM_H
#define G4DETECTORS_PHG4SPACALSUBSYSTEM_H

#include "PHG4DetectorSubsystem.h"

#include <string>                   // for string

class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4SpacalDetector;
class PHG4SteppingAction;

class PHG4SpacalSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4SpacalSubsystem(const std::string &name = "PHG4SpacalSubsystem",
                      const int layer = 0);

  //! destructor
  virtual ~PHG4SpacalSubsystem();

  //! init
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

  //! accessors (reimplemented)
  virtual PHG4Detector *GetDetector() const;
  virtual PHG4SteppingAction *GetSteppingAction() const { return steppingAction_; }

  PHG4DisplayAction *GetDisplayAction() const { return m_DisplayAction; }

  void
  Print(const std::string &what = "ALL") const;

 private:
  void SetDefaultParameters();
  //  SpacalGeom_t _geom;

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4SpacalDetector *detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction *steppingAction_;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction *m_DisplayAction;
};

#endif
