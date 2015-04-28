#ifndef PHG4VCylinderSteppingAction_h
#define PHG4VCylinderSteppingAction_h

#include "g4main/PHG4SteppingAction.h"
#include <string>

class PHG4CylinderDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4CylinderSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4CylinderSteppingAction( PHG4CylinderDetector* );

  //! destroctor
  virtual ~PHG4CylinderSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  void set_zmin(const float z) {zmin = z;}
  void set_zmax(const float z) {zmax = z;}

  private:

  //! pointer to the detector
  PHG4CylinderDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4Hit *hit;
  float zmin;
  float zmax;
};


#endif //__G4PHPHYTHIAREADER_H__
