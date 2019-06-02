// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.1 $$
 * \date $$Date: 2014/03/24 01:36:44 $$
 */

#ifndef G4DETECTORS_PHG4SPACALSTEPPINGACTION_H
#define G4DETECTORS_PHG4SPACALSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class PHCompositeNode;
class PHG4SpacalDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;

class PHG4SpacalSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  explicit PHG4SpacalSteppingAction(PHG4SpacalDetector *);

  //! destroctor
  virtual ~PHG4SpacalSteppingAction();

  //! stepping action
  virtual bool
  UserSteppingAction(const G4Step *, bool);

  //! reimplemented from base class
  virtual void
  SetInterfacePointers(PHCompositeNode *);

  double
  get_zmin();

  double
  get_zmax();

 private:
  //! pointer to the detector
  PHG4SpacalDetector *detector_;

  //! pointer to hit container
  PHG4HitContainer *hits_;
  PHG4HitContainer *absorberhits_;
  PHG4Hit *hit;
  PHG4HitContainer *savehitcontainer;
  PHG4Shower *saveshower;
  int savetrackid;
  int savepoststepstatus;
};

#endif  // PHG4VHcalSteppingAction_h
