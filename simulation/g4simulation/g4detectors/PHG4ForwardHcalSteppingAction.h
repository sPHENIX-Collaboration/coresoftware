#ifndef PHG4VForwardHcalSteppingAction_h
#define PHG4VForwardHcalSteppingAction_h

#include <g4main/PHG4SteppingAction.h>
#include <Geant4/G4Step.hh>


class PHG4ForwardHcalDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4ForwardHcalSteppingAction : public PHG4SteppingAction
{

public:

  //! constructor
  PHG4ForwardHcalSteppingAction( PHG4ForwardHcalDetector* );

  //! destroctor
  virtual ~PHG4ForwardHcalSteppingAction()
  {}

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
  PHG4HitContainer * hits_;
  PHG4HitContainer * absorberhits_;
  PHG4Hit *hit;

  int absorbertruth; 
  int light_scint_model; 

};


#endif // PHG4ForwardHcalSteppingAction_h
