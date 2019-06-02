// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PHENIXSTEPPINGACTION_H
#define G4MAIN_PHG4PHENIXSTEPPINGACTION_H

#include <Geant4/G4UserSteppingAction.hh>
#include <list>

class G4Step;
class PHG4SteppingAction;
class PHCompositeNode;

class PHG4PhenixSteppingAction : public G4UserSteppingAction
{

  public:
  PHG4PhenixSteppingAction( void )
  {}

  virtual ~PHG4PhenixSteppingAction();
  

  //! register an action. This is called in PHG4Reco::Init based on which actions are found on the tree
  void AddAction( PHG4SteppingAction* action )
  {
    if (action)
      {
	actions_.push_back( action );
      }
  }

  virtual void UserSteppingAction(const G4Step*);

  private:

  //! list of subsystem specific stepping actions
  typedef std::list<PHG4SteppingAction*> ActionList;
  ActionList actions_;

};


#endif
