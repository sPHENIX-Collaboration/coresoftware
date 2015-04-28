#ifndef PHG4VCylinderRegionSteppingAction_h
#define PHG4VCylinderRegionSteppingAction_h

#include <Geant4/G4UserSteppingAction.hh>


class PHCompositeNode;
class PHG4CylinderDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4CylinderRegionSteppingAction : public G4UserSteppingAction
{

  public:

  //! constructor
  PHG4CylinderRegionSteppingAction( PHG4CylinderDetector* );

  //! destroctor
  virtual ~PHG4CylinderRegionSteppingAction()
  {}

  //! stepping action
  virtual void UserSteppingAction(const G4Step*);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  private:

  //! pointer to the detector
  PHG4CylinderDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4Hit *hit;
};


#endif //__G4PHPHYTHIAREADER_H__
