// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4TRACKINGACTION_H
#define G4MAIN_PHG4TRACKINGACTION_H

#include <Geant4/G4UserTrackingAction.hh>

class G4Track;
class PHCompositeNode;

class PHG4TrackingAction : public G4UserTrackingAction
{
public:
  PHG4TrackingAction( void ) {}

  ~PHG4TrackingAction() override {}

//   //! tracking action. This defines pre/post processing of a single track in an event
  void PreUserTrackingAction(const G4Track*) override {}

   void PostUserTrackingAction(const G4Track*) override {}

  //! Set the node pointers
  virtual void SetInterfacePointers( PHCompositeNode* ) {return;}

  virtual int ResetEvent(PHCompositeNode *) {return 0;}

};


#endif // G4MAIN_PHG4TRACKINGACTION_H
