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
 * What the reconstruction reads from one event: the track map, the vertex maps the configuration asks for,
 * and the beam crossings of the intermediate candidates accepted so far.
 *
 * Built at the start of KFParticle_eventReconstruction::createDecay and passed down by reference; it lives for
 * that one call. KFParticle_sPHENIX::InitRun builds one to check that the nodes exist (missingNodes) and aborts
 * the run otherwise, so the map pointers are not null-checked here.
 */
class KFParticle_eventView
{
 public:
  /// Fetches the nodes the configuration needs: the track map; unless a fake PV is used, the MBD or SVTX
  /// vertex map and, unless dontUseGlobalVertex, the GlobalVertexMap
  KFParticle_eventView(PHCompositeNode *topNode, const KFParticle_Config &config);

  /// Nodes the configuration needs that were not found, as "<class> node <name>"; empty if all were found
  const std::vector<std::string> &missingNodes() const { return m_missingNodes; }

  PHCompositeNode *topNode() const { return m_topNode; }
  const SvtxTrackMap *trackMap() const { return m_trackMap; }
  const GlobalVertexMap *globalVertexMap() const { return m_globalVertexMap; }

  /// The MBD or SVTX vertex (whichever useMbdVertex selected) with this id; nullptr if it isn't there
  const Vertex *findDstVertex(unsigned int id) const;

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
