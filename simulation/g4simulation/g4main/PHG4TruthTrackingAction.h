#ifndef G4MAIN_PHG4TRUTHTRACKINGACTION_H
#define G4MAIN_PHG4TRUTHTRACKINGACTION_H

#include "PHG4TrackingAction.h"

#include <Geant4/G4ThreeVector.hh>

#include <map>

class PHG4TruthInfoContainer;
class PHG4TruthEventAction;

class PHG4TruthTrackingAction : public PHG4TrackingAction
{
 public:
  //! constructor
  PHG4TruthTrackingAction(PHG4TruthEventAction*);

  //! destructor
  virtual ~PHG4TruthTrackingAction() {}

  //! tracking action
  virtual void PreUserTrackingAction(const G4Track*);

  virtual void PostUserTrackingAction(const G4Track*);

  //! Set pointers to the i/o nodes
  virtual void SetInterfacePointers(PHCompositeNode*);

  int ResetEvent(PHCompositeNode*);

 private:
  std::map<G4ThreeVector, int> m_VertexMap;

  //! pointer to the "owning" event action
  PHG4TruthEventAction* m_EventAction;

  //! pointer to truth information container
  PHG4TruthInfoContainer* m_TruthInfoList;
};

#endif
