#include "TrackingIterationCounter.h"

#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>
//____________________________________________________________________________..
TrackingIterationCounter::TrackingIterationCounter(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
TrackingIterationCounter::~TrackingIterationCounter()
{
}

//____________________________________________________________________________..
int TrackingIterationCounter::Init(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TrackingIterationCounter::InitRun(PHCompositeNode *topNode)
{
  int val = getNodes(topNode);
  if (val != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  val = createNodes(topNode);
  return val;
}

//____________________________________________________________________________..
int TrackingIterationCounter::process_event(PHCompositeNode *)
{
  for (const auto &[key, track] : *m_trackMap)
  {
    auto silseed = track->get_silicon_seed();
    auto tpcseed = track->get_tpc_seed();

    addClustersToIterationMap(silseed);
    addClustersToIterationMap(tpcseed);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void TrackingIterationCounter::addClustersToIterationMap(TrackSeed *seed)
{
  for (auto clusIter = seed->begin_cluster_keys();
       clusIter != seed->end_cluster_keys();
       ++clusIter)
  {
    TrkrDefs::cluskey key = *clusIter;
    m_iterMap->addIteration(key, m_iteration);
  }
}

//____________________________________________________________________________..
int TrackingIterationCounter::End(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
int TrackingIterationCounter::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHActsTrkFitter::createNodes");
  }

  PHNodeIterator dstIter(topNode);
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_iterMap = findNode::getClass<TrkrClusterIterationMap>(topNode, "TrkrClusterIterationMap");
  if (!m_iterMap)
  {
    m_iterMap = new TrkrClusterIterationMapv1;
    auto node =
        new PHIODataNode<PHObject>(m_iterMap, "TrkrClusterIterationMap", "PHObject");
    svtxNode->addNode(node);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
int TrackingIterationCounter::getNodes(PHCompositeNode *topNode)
{
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!m_trackMap)
  {
    std::cout << PHWHERE << "No track map, bailing. " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
