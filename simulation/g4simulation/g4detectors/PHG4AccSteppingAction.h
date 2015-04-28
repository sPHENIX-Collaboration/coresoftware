#ifndef PHG4VAccSteppingAction_h
#define PHG4VAccSteppingAction_h

#include <g4main/PHG4SteppingAction.h>
#include <string>

class PHG4AccDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4AccSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4AccSteppingAction( PHG4AccDetector* );

  //! destroctor
  virtual ~PHG4AccSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  void set_zmin(const float z) {zmin = z;}
  void set_zmax(const float z) {zmax = z;}

  private:

  //! pointer to the detector
  PHG4AccDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4HitContainer * absorberhits_;
  PHG4Hit *hit;
  float zmin;
  float zmax;
};


#endif //__G4PHPHYTHIAREADER_H__
