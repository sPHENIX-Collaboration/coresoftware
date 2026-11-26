#ifndef KFPARTICLESPHENIX_KFPARTICLEEVENTDATA_H
#define KFPARTICLESPHENIX_KFPARTICLEEVENTDATA_H

#include "KFParticle_util.h"

#include <KFParticle.h>

#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

struct KFParticle_Config;
class GlobalVertexMap;
class PHCompositeNode;
class SvtxTrackMap;

/**
 * The DST nodes the reconstruction reads, fetched once per event (step 1 of 2).
 *
 * Which nodes are fetched depends on the configuration: always the track map; unless a fake PV is used, the MBD
 * or SVTX vertex map and, unless dontUseGlobalVertex, the GlobalVertexMap.
 *
 * KFParticle_sPHENIX::InitRun calls fetch() once just to check that the nodes exist (missing) and aborts the run
 * otherwise, so the node pointers are not checked again during the event.
 */
struct KFParticle_eventNodes
{
  const SvtxTrackMap *trackMap{nullptr};
  KFParticle_util::DstVertices dstVertices;         // the MBD or SVTX vertices (useMbdVertex) by id; empty with a fake PV
  const GlobalVertexMap *globalVertexMap{nullptr};  // nullptr with dontUseGlobalVertex or a fake PV
  std::vector<std::string> missing;                 // "<class> node <name>" for every node that was not found

  static KFParticle_eventNodes fetch(PHCompositeNode *topNode, const KFParticle_Config &config);
};

/**
 * Everything the candidate search needs to know about the current event, worked out once (step 2 of 2).
 *
 * Why it exists: many functions in the search (track selection, vertex matching, building mothers) need the same
 * facts about the event: which beam crossing a primary vertex or a track belongs to, whether a track is TPC-only,
 * a track's dE/dx. Instead of each function looking them up in the DST nodes again (many times per event),
 * createDecay works them out once into this struct and hands it down as a function argument. It holds plain
 * values only, no DST nodes.
 *
 * Lifetime: one per event. It is created at the top of KFParticle_eventReconstruction::createDecay and destroyed
 * when that function returns, so nothing carries over to the next event.
 *
 * Who may change it: functions that only read take `const KFParticle_eventData &`. Only buildMother writes
 * (registerIntermediate, when it accepts an intermediate candidate), so it and the functions that call it take a
 * non-const reference. (registerIntermediate is static so that this stays a plain struct with public fields.)
 *
 * Use (in createDecay):
 *   const KFParticle_eventNodes nodes = KFParticle_eventNodes::fetch(topNode, m_config);
 *   std::vector<KFParticle> daughterParticles = makeAllDaughterParticles(nodes.trackMap);
 *   KFParticle_eventData eventData = KFParticle_eventData::make(topNode, nodes, daughterParticles, m_config);
 */
struct KFParticle_eventData
{
  /// A primary vertex together with the beam crossing of the DST vertex at its position
  struct PrimaryVertex
  {
    KFParticle vertex;
    /// std::nullopt: no DST vertex was found for this PV, so no particle matches it.
    /// INT_MIN: DST vertices were found but none at this position; only particles without crossings match it.
    /// Not used with a fake PV, which matches every particle.
    std::optional<int> crossing;
  };

  /// Every primary vertex of the event (a single fake PV with useFakePrimaryVertex)
  std::vector<PrimaryVertex> primaryVertices;

  /// Beam crossings by particle id: every track in the track map has {its crossing}; every intermediate candidate
  /// accepted so far in this event has the crossings of its daughters, under -2, -3, ... (registerIntermediate)
  std::map<int, std::vector<int>> crossingsById;
  int nextIntermediateId{-2};

  /// Ids of the tracks in crossing 0 without silicon clusters (TPC-only)
  std::set<int> tpcOnlyAtZeroCrossing;

  /// usePID only: dE/dx of every daughter track (KFParticle_util::get_dEdx); not in here: not a track (-1)
  std::map<int, float> dEdxById;

  /// Builds all of the above from the fetched nodes; topNode is only used for the dE/dx
  static KFParticle_eventData make(PHCompositeNode *topNode, const KFParticle_eventNodes &nodes,
                                   const std::vector<KFParticle> &daughterParticles, const KFParticle_Config &config);

  /// Gives an accepted intermediate candidate the next ID (-2, -3, ... within the event) and records its beam crossings
  static int registerIntermediate(KFParticle_eventData &eventData, const std::vector<int> &crossings);
};

#endif  // KFPARTICLESPHENIX_KFPARTICLEEVENTDATA_H
