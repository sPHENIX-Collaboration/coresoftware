#ifndef KFPARTICLESPHENIX_KFPARTICLEEVENTVIEW_H
#define KFPARTICLESPHENIX_KFPARTICLEEVENTVIEW_H

#include "KFParticle_util.h"

#include <map>
#include <string>
#include <vector>

struct KFParticle_Config;
class GlobalVertexMap;
class PHCompositeNode;
class SvtxTrackMap;
class Vertex;

/**
 * Everything the candidate search needs from the current event, gathered in one object.
 *
 * Why it exists: many functions in the search (track selection, vertex matching, building mothers) need the
 * same few DST nodes. Instead of each function looking them up again, createDecay fetches them once into this
 * object and hands it down as a function argument.
 *
 * What it holds:
 *   1. DST nodes, fetched once in the constructor and only read afterwards: the track map, the MBD or SVTX
 *      vertices (copied into a map by id) and the GlobalVertexMap. Which of them are fetched depends on the
 *      configuration (see the constructor).
 *   2. One piece of bookkeeping that grows during the event: every intermediate candidate that is accepted gets
 *      an ID (-2, -3, ...) and its beam crossings are recorded under that ID, so that a later step can look up
 *      the crossings of a candidate built from it.
 *
 * Lifetime: one object per event. It is created at the top of KFParticle_eventReconstruction::createDecay and
 * destroyed when that function returns, so nothing carries over to the next event. Don't keep pointers into it.
 *
 * Who may change it: functions that only read take `const KFParticle_eventView &`. Only buildMother writes
 * (registerIntermediate), so it and the functions that call it take a non-const reference.
 *
 * Null checks: KFParticle_sPHENIX::InitRun builds one of these just to check that the nodes exist (missingNodes)
 * and aborts the run otherwise, so the node pointers are not checked again during the event.
 *
 * Use:
 *   KFParticle_eventView eventView(topNode, m_config);                         // fetch the nodes
 *   std::vector<KFParticle> tracks = makeAllDaughterParticles(eventView);       // read
 *   mother.SetId(eventView.registerIntermediate(crossings));                    // the only write
 */
class KFParticle_eventView
{
 public:
  /// Fetches the nodes the configuration needs: the track map; unless a fake PV is used, the MBD or SVTX
  /// vertex map and, unless dontUseGlobalVertex, the GlobalVertexMap
  KFParticle_eventView(PHCompositeNode *topNode, const KFParticle_Config &config);

  /// Nodes the configuration needs that were not found, as "<class> node <name>"; empty if all were found
  const std::vector<std::string> &missingNodes() const { return m_missingNodes; }

  // ---- 1. DST nodes (read only)

  PHCompositeNode *topNode() const { return m_topNode; }
  const SvtxTrackMap *trackMap() const { return m_trackMap; }
  const GlobalVertexMap *globalVertexMap() const { return m_globalVertexMap; }

  /// The MBD or SVTX vertices (whichever useMbdVertex selected) by id; empty with a fake PV
  const KFParticle_util::DstVertices &dstVertices() const { return m_dstVertices; }

  /// The MBD or SVTX vertex (whichever useMbdVertex selected) with this id; nullptr if it isn't there
  const Vertex *findDstVertex(unsigned int id) const;

  // ---- 2. Intermediate candidates accepted so far in this event

  /// Gives an accepted intermediate candidate the next ID (-2, -3, ... within the event) and records its beam crossings
  int registerIntermediate(const std::vector<int> &crossings);

  /// Beam crossings recorded for this intermediate ID; nullptr if no intermediate has this ID
  const std::vector<int> *intermediateCrossings(int id) const;

 private:
  PHCompositeNode *m_topNode{nullptr};

  SvtxTrackMap *m_trackMap{nullptr};
  KFParticle_util::DstVertices m_dstVertices;  // the MBD or the SVTX vertices, by id
  GlobalVertexMap *m_globalVertexMap{nullptr};
  std::vector<std::string> m_missingNodes;

  std::map<int, std::vector<int>> m_intermediateCrossings;
  int m_nextIntermediateId{-2};
};

#endif  // KFPARTICLESPHENIX_KFPARTICLEEVENTVIEW_H
