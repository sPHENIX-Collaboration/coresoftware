// $$Id: PHG4RICHSteppingAction.h,v 1.4 2013/12/22 19:33:38 nfeege Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov> and Nils Feege <nils.feege@stonybrook.edu>
 * \version $$Revision: 1.4 $$
 * \date $$Date: 2013/12/22 19:33:38 $$
 */

#ifndef PHG4RICHSteppingAction_h
#define PHG4RICHSteppingAction_h

#include <Geant4/G4UserSteppingAction.hh>

#include "Geant4/G4OpBoundaryProcess.hh"

class PHCompositeNode;
class PHG4RICHDetector;
class PHG4Hit;
class PHG4HitContainer;

/**
   * \brief This class defines the user stepping action for the ePHENIX RICH volumes
   * within Fun4All.
   *
   * Optical photons can be absorbed/detected at the optical surface of the photocathode.
   * The efficiency for detection is a parameter of this optical surface defined during
   * detector construction (set in ePHENIXRICH::RICH_Geometry). The x,y,z positions of
   * where the photons are detected are stored in a PHG4Hits collection.
   *
   * \see ePHENIXRICH::RICH_Geometry
   * \see ePHENIXRICH::ePHENIXRICHConstruction
   * \see PHG4RICHDetector
   * \see PHG4RICHSteppingAction
   * \see PHG4RICHSubsystem
   *
   */
class PHG4RICHSteppingAction : public G4UserSteppingAction
{
 public:
  PHG4RICHSteppingAction(PHG4RICHDetector*);
  virtual ~PHG4RICHSteppingAction() {}

  virtual void UserSteppingAction(const G4Step*);

  virtual void SetInterfacePointers(PHCompositeNode*);

 private:
  bool MakeHit(const G4Step* aStep);

  PHG4RICHDetector* detector_;
  PHG4HitContainer* hits_;
  PHG4Hit* hit;

  G4OpBoundaryProcessStatus fExpectedNextStatus;
};

#endif
