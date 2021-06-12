// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PHENIXSTACKINGACTION_H
#define G4MAIN_PHG4PHENIXSTACKINGACTION_H

#include <Geant4/G4UserStackingAction.hh>

#include <list>

class G4Step;
class PHG4StackingAction;

class PHG4PhenixStackingAction : public G4UserStackingAction
{

  public:
  PHG4PhenixStackingAction( void )
  {}

  ~PHG4PhenixStackingAction() override;
  

  //! register an action. This is called in PHG4Reco::Init based on which actions are found on the tree
  void AddAction( PHG4StackingAction* action )
  {
    if (action)
      {
	actions_.push_back( action );
      }
  }

  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack) override;
  void PrepareNewEvent() override;

  private:

  //! list of subsystem specific stacking actions
  typedef std::list<PHG4StackingAction*> ActionList;
  ActionList actions_;

};


#endif
