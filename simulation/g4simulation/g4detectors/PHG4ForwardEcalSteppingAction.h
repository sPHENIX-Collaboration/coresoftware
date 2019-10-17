// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4FORWARDECALSTEPPINGACTION_H
#define G4DETECTORS_PHG4FORWARDECALSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <Geant4/G4TouchableHandle.hh>

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4ForwardEcalDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;
class PHParameters;

class PHG4ForwardEcalSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4ForwardEcalSteppingAction(PHG4ForwardEcalDetector *, const PHParameters *parameters);

  //! destroctor
  virtual ~PHG4ForwardEcalSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*);

 private:
  int FindTowerIndex(G4TouchableHandle touch, int& j, int& k);

  int ParseG4VolumeName(G4VPhysicalVolume* volume, int& j, int& k);

  //! pointer to the detector
  PHG4ForwardEcalDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer* hits_;
  PHG4HitContainer* absorberhits_;
  const PHParameters *m_Params;
  PHG4HitContainer* hitcontainer;
  PHG4Hit* hit;
  PHG4Shower* saveshower;

  int absorbertruth;
  int light_scint_model;
  int m_IsBlackHole;
};

#endif  // PHG4ForwardEcalSteppingAction_h
