// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CONEREGIONSTEPPINGACTION_H
#define G4DETECTORS_PHG4CONEREGIONSTEPPINGACTION_H

#include <Geant4/G4UserSteppingAction.hh>


class G4Step;
class PHCompositeNode;
class PHG4ConeDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4ConeRegionSteppingAction : public G4UserSteppingAction
{

  public:

  //! constructor
  PHG4ConeRegionSteppingAction( PHG4ConeDetector* );

  //! destroctor
  virtual ~PHG4ConeRegionSteppingAction()
  {}

  //! stepping action
  virtual void UserSteppingAction(const G4Step*);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  private:

  //! pointer to the detector
  PHG4ConeDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4Hit *hit;
};


#endif //__G4PHPHYTHIAREADER_H__
