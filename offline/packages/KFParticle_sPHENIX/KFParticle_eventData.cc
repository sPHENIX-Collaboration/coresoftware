#include "KFParticle_eventData.h"

#include "KFParticle_Config.h"
#include "KFParticle_util.h"

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_v5.h>

#include <phool/getClass.h>

#include <cmath>   // for abs
#include <limits>

namespace
{
  // The node, or nullptr with "<type> node <name>" added to missing
  template <class T>
  T *fetchNode(PHCompositeNode *topNode, const std::string &type, const std::string &name, std::vector<std::string> &missing)
  {
    T *node = findNode::getClass<T>(topNode, name);
    if (!node)
    {
      missing.push_back(type + " node " + name);
    }
    return node;
  }

  // True if the two positions agree to within 1e-6 cm in every coordinate
  bool vertexToleranceCheck(float vertexMap_pos[3], float kfp_vertex_pos[3])
  {
    constexpr float epsilon = 1e-6F;
    for (unsigned int i = 0; i < 3; ++i)
    {
      if (std::abs(vertexMap_pos[i] - kfp_vertex_pos[i]) > epsilon)
      {
        return false;
      }
    }
    return true;
  }

  // One fake PV at the origin with no errors (useFakePrimaryVertex)
  KFParticle createFakePV()
  {
    float f_vertexParameters[6] = {0};

    float f_vertexCovariance[21] = {0};

    KFParticle kfp_vertex;
    kfp_vertex.Create(f_vertexParameters, f_vertexCovariance, 0, -1);
    kfp_vertex.NDF() = 0;
    kfp_vertex.Chi2() = 0;
    kfp_vertex.SetId(0);
    return kfp_vertex;
  }

  // The beam crossing of the DST vertex at the position of this PV (MBD: z only, x and y stay 0):
  //   direct (dontUseGlobalVertex): the DST vertex with the PV's id;
  //   global: every member of the selected type of the global vertex with the PV's id; the last one at the position wins.
  // std::nullopt if there is no DST vertex / global vertex / member of the type to look at; INT_MIN if none is at the position.
  std::optional<int> pvCrossing(const KFParticle &vertex, const KFParticle_eventNodes &nodes, const KFParticle_Config &config)
  {
    int vertexCrossing = std::numeric_limits<int>::min();

    float obtained_vertex[3] = {0, 0, std::numeric_limits<float>::lowest()};
    float kfp_vertex[3] = {vertex.GetX(), vertex.GetY(), vertex.GetZ()};

    auto matchCrossing = [&](const Vertex *dstVertex)
    {
      if (!KFParticle_util::isMbdVertex(dstVertex))
      {
        obtained_vertex[0] = dstVertex->get_x();
        obtained_vertex[1] = dstVertex->get_y();
      }
      obtained_vertex[2] = dstVertex->get_z();
      if (vertexToleranceCheck(obtained_vertex, kfp_vertex))
      {
        vertexCrossing = dstVertex->get_beam_crossing();
      }
    };

    if (config.dontUseGlobalVertex)
    {
      const Vertex *dstVertex = KFParticle_util::findDstVertex(nodes.dstVertices, vertex.Id());
      if (!dstVertex)
      {
        return std::nullopt;
      }
      matchCrossing(dstVertex);
    }
    else
    {
      const GlobalVertex *globalVertex = nodes.globalVertexMap->get(vertex.Id());
      if (!globalVertex)
      {
        return std::nullopt;
      }

      GlobalVertex::VTXTYPE whichVtx = config.useMbdVertex ? GlobalVertex::MBD : GlobalVertex::SVTX;

      auto svtxiter = globalVertex->find_vertexes(whichVtx);
      if (svtxiter == globalVertex->end_vertexes())
      {
        return std::nullopt;
      }
      const auto &svtxvertexvector = svtxiter->second;

      for (const auto &gvertex : svtxvertexvector)
      {
        const Vertex *dstVertex = KFParticle_util::findDstVertex(nodes.dstVertices, gvertex->get_id());
        if (!dstVertex)
        {
          continue;
        }
        matchCrossing(dstVertex);
      }
    }

    return vertexCrossing;
  }

  // In crossing 0 without silicon clusters
  bool isTPConlyAtZeroCrossing(const SvtxTrack *svtxTrack)
  {
    if (svtxTrack->get_crossing() != 0)
    {
      return false;
    }

    // SvtxTrack_v5 keeps no seed pointers: TPC-only means no silicon clusters on the track
    if (const auto *v5 = dynamic_cast<const SvtxTrack_v5 *>(svtxTrack))
    {
      return v5->get_nmvtx_clusters() == 0 && v5->get_nintt_clusters() == 0;
    }

    // SvtxTrack_v4 and earlier carry the seed pointer
    return svtxTrack->get_silicon_seed() == nullptr;
  }
}  // namespace

KFParticle_eventNodes KFParticle_eventNodes::fetch(PHCompositeNode *topNode, const KFParticle_Config &config)
{
  KFParticle_eventNodes nodes;

  nodes.trackMap = fetchNode<SvtxTrackMap>(topNode, "SvtxTrackMap", config.trackMapNodeName, nodes.missing);

  if (config.useFakePrimaryVertex)
  {
    return nodes;
  }

  if (config.useMbdVertex)
  {
    nodes.dstVertices = KFParticle_util::dstVertices(fetchNode<MbdVertexMap>(topNode, "MbdVertexMap", "MbdVertexMap", nodes.missing));
  }
  else
  {
    nodes.dstVertices = KFParticle_util::dstVertices(fetchNode<SvtxVertexMap>(topNode, "SvtxVertexMap", config.vertexMapNodeName, nodes.missing));
  }

  if (!config.dontUseGlobalVertex)
  {
    nodes.globalVertexMap = fetchNode<GlobalVertexMap>(topNode, "GlobalVertexMap", "GlobalVertexMap", nodes.missing);
  }

  return nodes;
}

KFParticle_eventData KFParticle_eventData::make(PHCompositeNode *topNode, const KFParticle_eventNodes &nodes,
                                                const std::vector<KFParticle> &daughterParticles, const KFParticle_Config &config)
{
  KFParticle_eventData eventData;

  if (config.useFakePrimaryVertex)
  {
    eventData.primaryVertices.push_back({createFakePV(), std::nullopt});
  }
  else
  {
    for (const KFParticle &vertex : KFParticle_util::makeAllPrimaryVertices(nodes.dstVertices, nodes.globalVertexMap, config))
    {
      eventData.primaryVertices.push_back({vertex, pvCrossing(vertex, nodes, config)});
    }
  }

  // every track, not only the daughters: the tables answer for any track id, as the track map did
  for (const auto &[id, track] : *nodes.trackMap)
  {
    eventData.crossingsById[static_cast<int>(id)] = {track->get_crossing()};
    if (isTPConlyAtZeroCrossing(track))
    {
      eventData.tpcOnlyAtZeroCrossing.insert(static_cast<int>(id));
    }
  }

  // only the daughters: buildMother asks for no other track, and get_dEdx is expensive
  if (config.usePID)
  {
    for (const KFParticle &daughter : daughterParticles)
    {
      eventData.dEdxById[daughter.Id()] = KFParticle_util::get_dEdx(topNode, daughter, config.trackMapNodeName);
    }
  }

  return eventData;
}

int KFParticle_eventData::registerIntermediate(KFParticle_eventData &eventData, const std::vector<int> &crossings)
{
  const int id = eventData.nextIntermediateId--;
  eventData.crossingsById[id] = crossings;
  return id;
}
