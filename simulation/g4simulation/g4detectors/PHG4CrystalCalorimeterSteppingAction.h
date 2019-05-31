#ifndef G4DETECTORS_PHG4CRYSTALCALORIMETERSTEPPINGACTION_H
#define G4DETECTORS_PHG4CRYSTALCALORIMETERSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>
#include <Geant4/G4Step.hh>

class PHG4CrystalCalorimeterDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;

class PHG4CrystalCalorimeterSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4CrystalCalorimeterSteppingAction(PHG4CrystalCalorimeterDetector*);

  //! destroctor
  virtual ~PHG4CrystalCalorimeterSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode*);

 private:
  //! Find tower index of mother volume
  int FindTowerIndex(G4TouchableHandle touch, int& j, int& k);

  //! Find combined tower index of volume, mother volume, and mother+1 volume
  int FindTowerIndex2LevelUp(G4TouchableHandle touch, int& j, int& k);

  int ParseG4VolumeName(G4VPhysicalVolume* volume, int& j, int& k);

  //! pointer to the detector
  PHG4CrystalCalorimeterDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer* hits_;
  PHG4HitContainer* absorberhits_;
  PHG4Hit* hit;
  PHG4HitContainer* savehitcontainer;
  PHG4Shower* saveshower;

  //  int light_scint_model;
};

#endif  // PHG4CrystalCalorimeterSteppingAction_h
