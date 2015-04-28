#ifndef PHG4VNullSteppingAction_h
#define PHG4VNullSteppingAction_h

#include "PHG4SteppingAction.h"

class PHG4Detector;

class PHG4NullSteppingAction : public PHG4SteppingAction
{
public:

  //! constructor
  PHG4NullSteppingAction( PHG4Detector* ) {};

  //! destroctor
  virtual ~PHG4NullSteppingAction()
  {}

  //! stepping action
  virtual bool UserSteppingAction(const G4Step*, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers( PHCompositeNode* ) {};

private:

};


#endif //__G4PHPHYTHIAREADER_H__
