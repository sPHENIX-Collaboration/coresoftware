// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PHENIXEVENTACTION_H
#define G4MAIN_PHG4PHENIXEVENTACTION_H

#include <phool/PHTimeServer.h>

#include <Geant4/G4UserEventAction.hh>

#include <list>


class G4Event;
class PHG4EventAction;

class PHG4PhenixEventAction : public G4UserEventAction
{

  public:
  PHG4PhenixEventAction( void );

  ~PHG4PhenixEventAction() override;

  //! register an action. This is called in PHG4Reco::Init based on which actions are found on the tree
  void AddAction( PHG4EventAction* action )
  { actions_.push_back( action ); }

  void BeginOfEventAction(const G4Event*) override;

  void EndOfEventAction(const G4Event*) override;

  private:

  //! list of subsystem specific Event actions
  typedef std::list<PHG4EventAction*> ActionList;
  ActionList actions_;

  //! module timer.
  PHTimeServer::timer _timer;
};


#endif
