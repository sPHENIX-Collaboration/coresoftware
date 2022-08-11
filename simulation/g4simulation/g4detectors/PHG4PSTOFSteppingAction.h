// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4PSTOFSTEPPINGACTION_H
#define G4DETECTORS_PHG4PSTOFSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4PSTOFDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHParametersContainer;

class PHG4PSTOFSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4PSTOFSteppingAction(PHG4PSTOFDetector*, const PHParametersContainer*);

  //! destructor
  ~PHG4PSTOFSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step*, bool) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode*) override;

 private:
  //! pointer to the detector
  PHG4PSTOFDetector* detector_ = nullptr;
  //! pointer to hit container
  PHG4HitContainer* hits_ = nullptr;
  PHG4Hit* hit = nullptr;
  PHG4HitContainer* savehitcontainer = nullptr;

  G4VPhysicalVolume* savevolpre = nullptr;
  G4VPhysicalVolume* savevolpost = nullptr;
  int savetrackid = -1;
  int saveprestepstatus = -1;
  int savepoststepstatus = -1;
  double edepsum = 0.;
  double eionsum = 0.;
};

#endif  // PHG4PSTOFSteppingAction_h__
