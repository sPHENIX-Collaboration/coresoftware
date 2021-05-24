// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PHENIXTRACKINGACTION_H
#define G4MAIN_PHG4PHENIXTRACKINGACTION_H

#include <Geant4/G4UserTrackingAction.hh>
#include <list>

// Master UserTrackingAction: tracking actions can be registered with this class
// and they will be called in the order of registration.

class G4Track;
class PHG4TrackingAction;

class PHG4PhenixTrackingAction : public G4UserTrackingAction
{
public:
  PHG4PhenixTrackingAction( void ) : verbosity_(0) {}

  ~PHG4PhenixTrackingAction() override;

  //! register an action. This is called in PHG4Reco::Init based on which actions are found on the tree
  void AddAction( PHG4TrackingAction* action ) { actions_.push_back( action ); }

  void PreUserTrackingAction(const G4Track*) override;

  void PostUserTrackingAction(const G4Track*) override;

  //! Get/Set verbosity level
  void Verbosity(int val) { verbosity_ = val; }
  int Verbosity() const { return verbosity_; }

private:

  //! list of subsystem specific Event actions
  typedef std::list<PHG4TrackingAction*> ActionList;
  ActionList actions_;
  int verbosity_;
};


#endif
