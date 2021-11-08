// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4TRUTHTRACKINGACTION_H
#define G4MAIN_PHG4TRUTHTRACKINGACTION_H

#include "PHG4TrackingAction.h"

#include <Geant4/G4ThreeVector.hh>

#include <map>
#include <vector>

class G4Track;
class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHG4TruthEventAction;
class PHG4Particle;
class PHG4VtxPoint;

class PHG4TruthTrackingAction : public PHG4TrackingAction
{
 public:
  //! constructor
  PHG4TruthTrackingAction(PHG4TruthEventAction*);

  //! destructor
  ~PHG4TruthTrackingAction() override {}

  //! tracking action
  void PreUserTrackingAction(const G4Track*) override;

  void PostUserTrackingAction(const G4Track*) override;

  //! Set pointers to the i/o nodes
  void SetInterfacePointers(PHCompositeNode*) override;

  int ResetEvent(PHCompositeNode*) override;

 private:
  std::map<G4ThreeVector, int> m_VertexMap;

  //! pointer to the "owning" event action
  PHG4TruthEventAction* m_EventAction;

  //! pointer to truth information container
  PHG4TruthInfoContainer* m_TruthInfoList;

  PHG4Particle* AddParticle(PHG4TruthInfoContainer&, G4Track&);
  PHG4VtxPoint* AddVertex(PHG4TruthInfoContainer&, const G4Track&);

  /// Machinery to keep track of upstream particles while adding Geant4 tracks
  /// to the truth info container
  ///@{
  void UpdateG4ParticleStack(const G4Track*);

  struct G4ParticleInfo { int g4track_id, particle_id, vertex_id; };
  std::vector<G4ParticleInfo> m_G4ParticleStack;
  G4ParticleInfo m_CurrG4Particle;
  ///@}
};

#endif
