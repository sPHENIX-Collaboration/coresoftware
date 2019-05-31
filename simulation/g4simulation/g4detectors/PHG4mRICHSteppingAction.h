#ifndef G4DETECTORS_PHG4MRICHSTEPPINGACTION_H
#define G4DETECTORS_PHG4MRICHSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class PHG4mRICHDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHParameters;
class PHG4Shower;
class G4VPhysicalVolume;

class PHG4mRICHSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4mRICHSteppingAction(PHG4mRICHDetector* detector,PHParameters* params );

  //! destroctor
  virtual ~PHG4mRICHSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  private:

  //! pointer to the detector
  PHG4mRICHDetector* detector_;

  //detector parameters
  // int active;
  int IsBlackHole;
  // int use_g4_steps;
  std::string detectorname;
  std::string superdetector;

  //! pointer to hit container
  PHG4HitContainer* hits_;
  PHG4HitContainer* absorberhits_;
  PHG4Hit* hit;
  PHG4HitContainer *savehitcontainer;
  PHG4Shower *saveshower;
  int savetrackid;
  int savepoststepstatus;

  int GetModuleID(G4VPhysicalVolume* volume);
};


#endif // PHG4mRICHSteppingAction_h
