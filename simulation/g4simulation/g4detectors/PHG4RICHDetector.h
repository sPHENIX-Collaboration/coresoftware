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

#include <string>                     // for string

class G4LogicalVolume;
class G4UserSteppingAction;
class PHCompositeNode;
class PHG4Subsystem;

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
  PHG4RICHDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const std::string& dname, const ePHENIXRICH::RICH_Geometry& g);
  PHG4RICHDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const std::string& dname);

  virtual ~PHG4RICHDetector(void)
  {
  }

  virtual void
  ConstructMe(G4LogicalVolume* world);

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
  G4Region* _region;
};

#endif
