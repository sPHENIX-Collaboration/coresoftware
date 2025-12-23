#include "TrackingIterationCounter.h"

#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/TrackSeedContainer.h>
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
int TrackingIterationCounter::process_event(PHCompositeNode *topNode)
{
  if (m_iterateSeeds)
  {
    iterateSeeds(topNode);

    return Fun4AllReturnCodes::EVENT_OK;
  }


  for (const auto &[key, track] : *m_trackMap)
  {
    auto *silseed = track->get_silicon_seed();
    auto *tpcseed = track->get_tpc_seed();
    if (silseed)
    {
     
        addClustersToIterationMap(silseed);
      
    }
    if (tpcseed)
    {
      addClustersToIterationMap(tpcseed);
    }
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
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  if (!m_trackMap)
  {
    if(!m_iterateSeeds){
    std::cout << PHWHERE << "No track map, bailing. " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void TrackingIterationCounter::iterateSeeds(PHCompositeNode *topNode)
{
  if (m_trackMapName.find("Silicon") != std::string::npos)
  {
    iterateSiliconSeeds(topNode);
  }
}
void TrackingIterationCounter::iterateSiliconSeeds(PHCompositeNode *topNode)
{
  auto seeds = findNode::getClass<TrackSeedContainer>(topNode, m_trackMapName);
  for (const auto &seed : *seeds)
  {
    if (!seed)
    {
      continue;
    }

    if (m_iteration == 1)
    {
      if (seed->size_cluster_keys() > 4)
      {
        int nmaps = 0;
        int nintt = 0;
        for (auto clusIter = seed->begin_cluster_keys();
             clusIter != seed->end_cluster_keys();
             ++clusIter)
        {
          auto trkrid = TrkrDefs::getTrkrId(*clusIter);
          if (trkrid == TrkrDefs::inttId)
          {
            nintt++;
          }
          if (trkrid == TrkrDefs::mvtxId)
          {
            nmaps++;
          }
        }
        if (nintt > 1 && nmaps > 2)
        {
          addClustersToIterationMap(seed);
        }
      }
    }
    else{
      addClustersToIterationMap(seed);
    }
  }
}
