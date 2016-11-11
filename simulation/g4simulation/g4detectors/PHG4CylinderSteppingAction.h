#ifndef PHG4CylinderSteppingAction_h
#define PHG4CylinderSteppingAction_h

#include <g4main/PHG4SteppingAction.h>

#include <string>

class PHG4CylinderDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;
class PHG4Parameters;

class PHG4CylinderSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4CylinderSteppingAction( PHG4CylinderDetector*, const PHG4Parameters *parameters );

  //! destroctor
  virtual ~PHG4CylinderSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  private:

  //! pointer to the detector
  PHG4CylinderDetector* detector_;

  const PHG4Parameters *params;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4Hit *hit;
  PHG4Shower *saveshower;
  int active;
  int IsBlackHole;
  double zmin;
  double zmax;
  double tmin;
  double tmax;
};

#endif
