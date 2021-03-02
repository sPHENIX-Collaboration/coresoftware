// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4FORWARDHCALSTEPPINGACTION_H
#define G4DETECTORS_PHG4FORWARDHCALSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <Geant4/G4TouchableHandle.hh>

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4ForwardHcalDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;
class PHParameters;

class PHG4ForwardHcalSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4ForwardHcalSteppingAction(PHG4ForwardHcalDetector*, const PHParameters* parameters);

  //! destructor
  virtual ~PHG4ForwardHcalSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*);

 private:
  int FindTowerIndex(G4TouchableHandle touch, int& j, int& k);

  int ParseG4VolumeName(G4VPhysicalVolume* volume, int& j, int& k);

  //! pointer to the detector
  PHG4ForwardHcalDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer* hits_;
  PHG4HitContainer* absorberhits_;
  PHG4HitContainer* hitcontainer;
  PHG4Hit* hit;
  PHG4Shower* saveshower;

  int m_IsActiveFlag = 0;
  int absorbertruth = 0;
  int light_scint_model;
  int m_IsBlackHole = 0;
};

#endif  // G4DETECTORS_PHG4FORWARDHCALSTEPPINGACTION_H
