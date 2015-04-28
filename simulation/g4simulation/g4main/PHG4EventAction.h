#ifndef PHG4EventAction_h
#define PHG4EventAction_h

#include <phool/PHCompositeNode.h>

class G4Event;

class PHG4EventAction
{

  public:
  PHG4EventAction( void )
  {}

  virtual ~PHG4EventAction()
  {}

  virtual void BeginOfEventAction(const G4Event*) {}

  virtual void EndOfEventAction(const G4Event*) {}

  //! get relevant nodes from top node passed as argument
  virtual void SetInterfacePointers( PHCompositeNode* ) {}

  virtual int ResetEvent(PHCompositeNode *) {return 0;}

};


#endif
