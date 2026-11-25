#include "KFParticle_eventView.h"

#include "KFParticle_Config.h"
#include "KFParticle_util.h"

#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <phool/getClass.h>

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
}  // namespace

KFParticle_eventView::KFParticle_eventView(PHCompositeNode *topNode, const KFParticle_Config &config)
  : m_topNode(topNode)
{
  m_trackMap = fetchNode<SvtxTrackMap>(topNode, "SvtxTrackMap", config.trackMapNodeName, m_missingNodes);

  if (config.useFakePrimaryVertex)
  {
    return;
  }

  if (config.useMbdVertex)
  {
    m_dstVertices = KFParticle_util::dstVertices(fetchNode<MbdVertexMap>(topNode, "MbdVertexMap", "MbdVertexMap", m_missingNodes));
  }
  else
  {
    m_dstVertices = KFParticle_util::dstVertices(fetchNode<SvtxVertexMap>(topNode, "SvtxVertexMap", config.vertexMapNodeName, m_missingNodes));
  }

  if (!config.dontUseGlobalVertex)
  {
    m_globalVertexMap = fetchNode<GlobalVertexMap>(topNode, "GlobalVertexMap", "GlobalVertexMap", m_missingNodes);
  }
}

const Vertex *KFParticle_eventView::findDstVertex(unsigned int id) const
{
  return KFParticle_util::findDstVertex(m_dstVertices, id);
}

int KFParticle_eventView::registerIntermediate(const std::vector<int> &crossings)
{
  const int id = m_nextIntermediateId--;
  m_intermediateCrossings[id] = crossings;
  return id;
}

const std::vector<int> *KFParticle_eventView::intermediateCrossings(int id) const
{
  auto it = m_intermediateCrossings.find(id);
  return it == m_intermediateCrossings.end() ? nullptr : &it->second;
}
