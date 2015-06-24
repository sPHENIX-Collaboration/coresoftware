#ifndef PHG4VCrystalCalorimeterSteppingAction_h
#define PHG4VCrystalCalorimeterSteppingAction_h

#include <g4main/PHG4SteppingAction.h>
#include <Geant4/G4Step.hh>


class PHG4CrystalCalorimeterDetector;
class PHG4Hitv8;
class PHG4HitContainer;

class PHG4CrystalCalorimeterSteppingAction : public PHG4SteppingAction
{

public:

  //! constructor
  PHG4CrystalCalorimeterSteppingAction( PHG4CrystalCalorimeterDetector* );

  //! destroctor
  virtual ~PHG4CrystalCalorimeterSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  int WhatAreYou(G4TouchableHandle touch, int& j, int& k);

  int ParseName(G4VPhysicalVolume* volume, int& j, int& k);

private:

  //! pointer to the detector
  PHG4CrystalCalorimeterDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4HitContainer * absorberhits_;
  PHG4Hitv8 *hit;

};


#endif // PHG4CrystalCalorimeterSteppingAction_h
