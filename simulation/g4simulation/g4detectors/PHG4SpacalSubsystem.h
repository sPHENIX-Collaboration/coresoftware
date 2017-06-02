/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.2 $$
 * \date $$Date: 2014/08/12 03:49:12 $$
 */

#ifndef PHG4SpacalSubsystem_h
#define PHG4SpacalSubsystem_h

#include "PHG4DetectorSubsystem.h"

#include <Geant4/G4String.hh>
#include <Geant4/G4Types.hh>

class PHG4SpacalDetector;
class PHG4SpacalSteppingAction;

class PHG4SpacalSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4SpacalSubsystem(const std::string &name = "PHG4SpacalSubsystem",
                      const int layer = 0);

  //! destructor
  virtual ~PHG4SpacalSubsystem(void)
  {
  }

  //! init
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

  //! accessors (reimplemented)
  virtual PHG4Detector *
  GetDetector(void) const;
  virtual PHG4SteppingAction *
  GetSteppingAction(void) const;

  void
  Print(const std::string &what = "ALL") const;

 private:
  //  SpacalGeom_t _geom;

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4SpacalDetector *detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SpacalSteppingAction *steppingAction_;

  void SetDefaultParameters();
};

#endif
