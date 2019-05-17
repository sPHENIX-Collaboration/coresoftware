// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4VFORWARDHCALSTEPPINGACTION_H
#define G4DETECTORS_PHG4VFORWARDHCALSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <Geant4/G4TouchableHandle.hh>

class G4Step;
class PHG4ForwardHcalDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;

class PHG4ForwardHcalSteppingAction : public PHG4SteppingAction
{

public:

  //! constructor
  PHG4ForwardHcalSteppingAction( PHG4ForwardHcalDetector* );

  //! destructor
  virtual ~PHG4ForwardHcalSteppingAction();


  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

private:

  int FindTowerIndex(G4TouchableHandle touch, int& j, int& k);

  int ParseG4VolumeName(G4VPhysicalVolume* volume, int& j, int& k);

  //! pointer to the detector
  PHG4ForwardHcalDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer *hits_;
  PHG4HitContainer *absorberhits_;
  PHG4HitContainer *hitcontainer;
  PHG4Hit *hit;
  PHG4Shower *saveshower;

  int absorbertruth; 
  int light_scint_model; 

};


#endif // G4DETECTORS_PHG4FORWARDHCALSTEPPINGACTION_H
