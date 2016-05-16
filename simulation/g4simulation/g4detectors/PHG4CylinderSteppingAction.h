#ifndef PHG4CylinderSteppingAction_h
#define PHG4CylinderSteppingAction_h

#include <g4main/PHG4SteppingAction.h>

#include <string>

class PHG4CylinderDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;

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

  void flush_cached_values();

  private:

  void save_previous_g4hit();

  //! pointer to the detector
  PHG4CylinderDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4Hit *hit;
  PHG4HitContainer *savehitcontainer;
  PHG4Shower *saveshower;
  int save_layer_id;
  float zmin;
  float zmax;
};


#endif
