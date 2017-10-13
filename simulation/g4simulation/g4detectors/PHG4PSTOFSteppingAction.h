#ifndef PHG4PSTOFSteppingAction_h__
#define PHG4PSTOFSteppingAction_h__

#include <g4main/PHG4SteppingAction.h>

class G4VPhysicalVolume;
class PHG4PSTOFDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4ParametersContainer;

class PHG4PSTOFSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4PSTOFSteppingAction(PHG4PSTOFDetector*, const PHG4ParametersContainer*);

  //! destructor
  virtual ~PHG4PSTOFSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*);

 private:
  //! pointer to the detector
  PHG4PSTOFDetector* detector_;
  const PHG4ParametersContainer* paramscontainer;
  //! pointer to hit container
  PHG4HitContainer* hits_;
  PHG4Hit* hit;
  PHG4HitContainer *savehitcontainer;

  G4VPhysicalVolume *savevolpre;
  G4VPhysicalVolume *savevolpost;
  int savetrackid;
  int saveprestepstatus;
  int savepoststepstatus;
  double edepsum;
  double eionsum;
};

#endif  // PHG4PSTOFSteppingAction_h__
