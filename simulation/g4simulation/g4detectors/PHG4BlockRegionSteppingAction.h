#ifndef PHG4VBlockRegionSteppingAction_h
#define PHG4VBlockRegionSteppingAction_h

#include <Geant4/G4UserSteppingAction.hh>


class PHCompositeNode;
class PHG4BlockDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4BlockRegionSteppingAction : public G4UserSteppingAction
{

  public:

  //! constructor
  PHG4BlockRegionSteppingAction( PHG4BlockDetector* );

  //! destroctor
  virtual ~PHG4BlockRegionSteppingAction()
  {}

  //! stepping action
  virtual void UserSteppingAction(const G4Step*);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  private:

  //! pointer to the detector
  PHG4BlockDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4Hit *hit;
};


#endif //__G4PHPHYTHIAREADER_H__
