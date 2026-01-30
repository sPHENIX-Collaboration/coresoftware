// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4TRUTHTRACKINGACTION_H
#define G4MAIN_PHG4TRUTHTRACKINGACTION_H

#include "PHG4TrackingAction.h"
#include "PHG4MCProcessDefs.h"

#include <Geant4/G4ThreeVector.hh>

#include <map>
#include <utility>
#include <vector>

class G4Track;
class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHG4TruthEventAction;
class PHG4Particle;
class PHG4VtxPoint;

/**
 * Construct a PHG4TruthTrackingAction associated with an event action.
 * @param eventAction Pointer to the owning PHG4TruthEventAction used to record per-event truth information.
 */

/**
 * Destroy the PHG4TruthTrackingAction.
 */

/**
 * Handle actions to perform before Geant4 begins tracking a G4 track.
 * @param track The Geant4 track about to be processed.
 */

/**
 * Handle actions to perform after Geant4 finishes tracking a G4 track.
 * @param track The Geant4 track that has just been processed.
 */

/**
 * Set required node/interface pointers from the given top-level node.
 * @param topNode Pointer to the PHCompositeNode root from which required I/O nodes are retrieved.
 * @returns Zero on success, non-zero on failure.
 */

/**
 * Reset per-event state using nodes found under the given composite node.
 * @param topNode Pointer to the PHCompositeNode for the current event.
 * @returns Zero on success, non-zero on failure.
 */

/**
 * Create or update a truth particle entry corresponding to the provided Geant4 track.
 * @param truth Container to which the particle entry will be added or updated.
 * @param track Geant4 track from which particle information is derived.
 * @returns Pointer to the created or updated PHG4Particle.
 */

/**
 * Create or update a truth vertex entry corresponding to the provided Geant4 track.
 * @param truth Container to which the vertex entry will be added or updated.
 * @param track Geant4 track whose production point will be recorded as a vertex.
 * @returns Pointer to the created or updated PHG4VtxPoint.
 */

/**
 * Determine whether a particle type is considered long-lived for truth-building.
 * @param pid Particle PDG identifier.
 * @returns `true` if the particle with the given PDG id is treated as long-lived, `false` otherwise.
 */

/**
 * Determine whether a particle should be flagged as an sPHENIX primary.
 * @param truth Truth information container used to evaluate primary status.
 * @param particle Particle to evaluate.
 * @returns `true` if the particle is considered an sPHENIX primary, `false` otherwise.
 */

/**
 * Update the internal upstream G4 particle stack when processing a new Geant4 track.
 * @param track Geant4 track used to update parent/ancestor particle bookkeeping.
 */
class PHG4TruthTrackingAction : public PHG4TrackingAction
{
 public:
  //! constructor
  PHG4TruthTrackingAction(PHG4TruthEventAction*);

  //! destructor
  ~PHG4TruthTrackingAction() = default;

  //! tracking action
  void PreUserTrackingAction(const G4Track*) override;

  void PostUserTrackingAction(const G4Track*) override;

  //! Set pointers to the i/o nodes
  void SetInterfacePointers(PHCompositeNode*) override;

  int ResetEvent(PHCompositeNode*) override;

 private:
  // Key is (position, process) to distinguish vertices at the same location but different processes
  std::map<std::pair<G4ThreeVector, PHG4MCProcess>, int> m_VertexMap;

  //! pointer to the "owning" event action
  PHG4TruthEventAction* m_EventAction;

  //! pointer to truth information container
  PHG4TruthInfoContainer* m_TruthInfoList{nullptr};

  PHG4Particle* AddParticle(PHG4TruthInfoContainer&, G4Track&);
  PHG4VtxPoint* AddVertex(PHG4TruthInfoContainer&, const G4Track&);

  // check if track is long-lived
  bool isLongLived(int pid) const;

  // check if track is sPHENIX primary
  bool issPHENIXPrimary(PHG4TruthInfoContainer& truth, PHG4Particle* particle) const;

  /// Machinery to keep track of upstream particles while adding Geant4 tracks
  /// to the truth info container
  ///@{
  void UpdateG4ParticleStack(const G4Track*);

  struct G4ParticleInfo
  {
    int g4track_id{0};
    int particle_id{0};
    int vertex_id{0};
  };
  std::vector<G4ParticleInfo> m_G4ParticleStack;
  G4ParticleInfo m_CurrG4Particle;
  ///@}
};

#endif
