// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4EVENTACTION_H
#define G4MAIN_PHG4EVENTACTION_H

class G4Event;
class PHCompositeNode;

class PHG4EventAction
{
 public:
  PHG4EventAction() = default;

  virtual ~PHG4EventAction() = default;

  virtual void BeginOfEventAction(const G4Event *) { return; }

  virtual void EndOfEventAction(const G4Event *) { return; }

  //! get relevant nodes from top node passed as argument
  virtual void SetInterfacePointers(PHCompositeNode *) { return; }

  virtual int ResetEvent(PHCompositeNode *) { return 0; }
};

#endif  // G4MAIN_PHG4EVENTACTION_H
