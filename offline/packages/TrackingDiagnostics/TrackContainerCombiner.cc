
#include "TrackContainerCombiner.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

//____________________________________________________________________________..
TrackContainerCombiner::TrackContainerCombiner(const std::string &name):
 SubsysReco(name)
{
}

//____________________________________________________________________________..
TrackContainerCombiner::~TrackContainerCombiner()
{
  
}

//____________________________________________________________________________..
int TrackContainerCombiner::InitRun(PHCompositeNode *topNode)
{
  int ret = getNodes(topNode);
  return ret;
}

//____________________________________________________________________________..
int TrackContainerCombiner::process_event(PHCompositeNode*)
{
  if(m_seedContainer)
  {
    if(Verbosity() > 1)
    {
      std::cout << "Seed container size to start " << m_newSeedContainer->size()
                << std::endl;
    }
    mergeSeeds();
    if (Verbosity() > 1)
    {
      std::cout << "Seed container size to finish " << m_newSeedContainer->size()
                << std::endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void TrackContainerCombiner::mergeSeeds()
{
  for (const auto& seed : *m_oldSeedContainer)
  {
      if(!seed)
      {
        continue;
      }
    m_newSeedContainer->insert(seed);
  }
}
//____________________________________________________________________________..
int TrackContainerCombiner::End(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrackContainerCombiner::getNodes(PHCompositeNode* topNode)
{
  m_newSeedContainer = findNode::getClass<TrackSeedContainer>(topNode, m_newContainerName);
  if(!m_newSeedContainer)
  {
    std::cout << PHWHERE << "No new track seed container, bailing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  m_oldSeedContainer = findNode::getClass<TrackSeedContainer>(topNode, m_oldContainerName);
  if (!m_oldSeedContainer)
  {
    std::cout << PHWHERE << "No old track seed container, bailing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}