// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.2 $$
 * \date $$Date: 2013/12/22 19:33:38 $$
 */

#ifndef G4DETECTORS_PHG4RICHDETECTOR_H
#define G4DETECTORS_PHG4RICHDETECTOR_H

#include "ePHENIXRICHConstruction.h"

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Region.hh>
#include <Geant4/G4Types.hh>
#include <Geant4/globals.hh>

#include <map>

class G4Material;
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class PHG4RICHDisplayAction;
class PHG4RICHSubsystem;

/**
   * \brief This class creates the ePHENIX RICH volumes for Geant4 within Fun4All via
   * ePHENIXRICH::ePHENIXRICHConstruction based on the geometry information in 
   * ePHENIXRICH::RICH_Geometry.
   *
   * \see ePHENIXRICH::RICH_Geometry
   * \see ePHENIXRICH::ePHENIXRICHConstruction
   * \see PHG4RICHDetector
   * \see PHG4RICHSteppingAction
   * \see PHG4RICHSubsystem
   *
   */
class PHG4RICHDetector : public PHG4Detector,
                         public ePHENIXRICH::ePHENIXRICHConstruction
{
 public:
  PHG4RICHDetector(PHG4RICHSubsystem* subsys, PHCompositeNode* Node, const ePHENIXRICH::RICH_Geometry& g);
  PHG4RICHDetector(PHG4RICHSubsystem* subsys, PHCompositeNode* Node);

  virtual ~PHG4RICHDetector(void)
  {
  }

  virtual void
  Construct(G4LogicalVolume* world);

  virtual G4UserSteppingAction*
  GetSteppingAction()
  {
    if (_region)
      return _region->GetRegionalSteppingAction();
    else
      return 0;
  }

  virtual void OverlapCheck(const bool chk = true)
  {
    PHG4Detector::OverlapCheck(chk);
    ePHENIXRICHConstruction::OverlapCheck(chk);
  }

 private:
  G4UserSteppingAction* stepping_action;

  G4Region* _region;
};

#endif
