// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4EVENTACTION_H
#define G4MAIN_PHG4EVENTACTION_H

class G4Event;
class PHCompositeNode;

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


#endif // G4MAIN_PHG4EVENTACTION_H
