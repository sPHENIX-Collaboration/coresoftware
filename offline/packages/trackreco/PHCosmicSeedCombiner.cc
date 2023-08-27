
#include "PHCosmicSeedCombiner.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase_historic/CosmicTrackSeed_v1.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>

#include <phool/PHCompositeNode.h>

//____________________________________________________________________________..
PHCosmicSeedCombiner::PHCosmicSeedCombiner(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
PHCosmicSeedCombiner::~PHCosmicSeedCombiner()
{
}

//____________________________________________________________________________..
int PHCosmicSeedCombiner::Init(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHCosmicSeedCombiner::InitRun(PHCompositeNode* topNode)
{
  int ret = getNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  return createNodes(topNode);
}

//____________________________________________________________________________..
int PHCosmicSeedCombiner::process_event(PHCompositeNode*)
{
  for (auto trackiter = m_seedMap->begin(); trackiter != m_seedMap->end();
       ++trackiter)
  {
    TrackSeed* track1 = *trackiter;
    if (!track1)
    {
      continue;
    }

    unsigned int tpcid1 = track1->get_tpc_seed_index();
    unsigned int siid1 = track1->get_silicon_seed_index();

    if (Verbosity() > 1)
    {
      std::cout << "tpc id" << tpcid1 << std::endl;
    }

    auto tpcseed1 = m_tpcSeeds->get(tpcid1);
    const float phi1 = tpcseed1->get_phi(m_clusterContainer, m_tGeometry);
    const float eta1 = tpcseed1->get_eta();

    for (auto trackiter2 = trackiter; trackiter2 != m_seedMap->end(); ++trackiter2)
    {
      if (trackiter == trackiter2) continue;
      TrackSeed* track2 = *trackiter2;
      if (!track2)
      {
        continue;
      }

      unsigned int tpcid2 = track2->get_tpc_seed_index();
      unsigned int siid2 = track2->get_silicon_seed_index();

      if (Verbosity() > 1)
      {
        std::cout << "tpc 2 " << tpcid2 << std::endl;
      }

      auto tpcseed2 = m_tpcSeeds->get(tpcid2);

      const float phi2 = tpcseed2->get_phi(m_clusterContainer, m_tGeometry);
      const float eta2 = tpcseed2->get_eta();

      //! phis are the same, so we subtract to check
      float dphi = phi1 - phi2;
      if (dphi > M_PI)
      {
        dphi -= 2. * M_PI;
      }
      else if (dphi < -1 * M_PI)
      {
        dphi += 2. * M_PI;
      }
      //! If they are back to back, dphi=pi
      dphi = fabs(dphi) - M_PI;

      //! etas are opposite each other, so we add them
      const float deta = eta1 + eta2;

      if (fabs(dphi) < 0.02 && fabs(deta) < 0.01)
      {
        auto seed = std::make_unique<CosmicTrackSeed_v1>();
        seed->set_silicon_seed_index(siid1);
        seed->set_silicon_seed_index2(siid2);
        seed->set_tpc_seed_index(tpcid1);
        seed->set_tpc_seed_index2(tpcid2);

        m_cosmicContainer->insert(seed.get());
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHCosmicSeedCombiner::End(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
int PHCosmicSeedCombiner::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHSiliconHelicalPropagator::createNodes");
  }

  PHNodeIterator dstIter(dstNode);
  PHCompositeNode* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_cosmicContainer = findNode::getClass<TrackSeedContainer>(topNode, "CosmicSeedContainer");
  if (!m_cosmicContainer)
  {
    m_cosmicContainer = new TrackSeedContainer_v1;
    PHIODataNode<PHObject>* trackNode = new PHIODataNode<PHObject>(m_cosmicContainer, "CosmicSeedContainer", "PHObject");
    svtxNode->addNode(trackNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
int PHCosmicSeedCombiner::getNodes(PHCompositeNode* topNode)
{
  m_tpcSeeds = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!m_tpcSeeds)
  {
    std::cout << PHWHERE << "TpcTrackSeedContainer not on node tree. Bailing"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_siliconSeeds = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if (!m_siliconSeeds)
  {
    std::cout << PHWHERE << "SiliconTrackSeedContainer not on node tree. Bailing"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterContainer)
  {
    std::cout << PHWHERE
              << "No trkr cluster container, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << "ActsGeometry not on node tree. Exiting."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_seedMap = findNode::getClass<TrackSeedContainer>(topNode, "SvtxTrackSeedContainer");
  if (!m_seedMap)
  {
    std::cout << "No Svtx seed map on node tree. Exiting."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
