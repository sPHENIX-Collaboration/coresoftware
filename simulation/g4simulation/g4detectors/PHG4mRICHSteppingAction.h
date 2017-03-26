#ifndef PHG4mRICHSteppingAction_h
#define PHG4mRICHSteppingAction_h

#include "g4main/PHG4SteppingAction.h"
//#include "PHG4Parameters.h"

class PHG4mRICHDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Parameters;

class PHG4mRICHSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4mRICHSteppingAction(PHG4mRICHDetector* detector,PHG4Parameters* params );

  //! destroctor
  virtual ~PHG4mRICHSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* );

  private:

  //! pointer to the detector
  PHG4mRICHDetector* detector_;

  //detector parameters
  int active;
  int IsBlackHole;
  int use_g4_steps;
  std::string detectorname;
  std::string superdetector;

  //! pointer to hit container
  PHG4HitContainer * hits_;
  PHG4HitContainer * absorberhits_;
  PHG4Hit *hit;
};


#endif // PHG4mRICHSteppingAction_h
