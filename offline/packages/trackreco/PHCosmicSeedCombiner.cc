
#include "PHCosmicSeedCombiner.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase_historic/TrackSeed.h>
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
  return getNodes(topNode);
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
    auto silseed1 = m_siliconSeeds->get(siid1);

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
      auto silseed2 = m_siliconSeeds->get(siid2);
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
      if (Verbosity() > 3)
      {
        std::cout << "phi 1 and phi2  " << phi1 << " , " << phi2 << std::endl;
        std::cout << "eta 1 and eta2 " << eta1 << " , " << eta2 << std::endl;
        std::cout << "dphi and deta " << dphi << " , " << deta << std::endl;
      }
      if (fabs(dphi) < m_dphiCut && fabs(deta) < m_detaCut)
      {
        //! add the clusters to the tpc seed and delete seed 2 since it is
        //! from the same track
        addKeys(tpcseed1, tpcseed2);
        if (silseed1)
        {
          if (silseed2)
          {
            addKeys(silseed1, silseed2);
          }
        }
        else
        {
          if (silseed2)
          {
            track1->set_silicon_seed_index(siid2);
          }
        }

        m_seedMap->erase(m_seedMap->index(trackiter2));
      }
      if (Verbosity() > 3)
      {
        track1->identify();
        tpcseed1->identify();
        if (silseed1)
        {
          silseed1->identify();
        }
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
void PHCosmicSeedCombiner::addKeys(TrackSeed* seedToAddTo, TrackSeed* seedToAdd)
{
  for (auto citer = seedToAdd->begin_cluster_keys();
       citer != seedToAdd->end_cluster_keys();
       ++citer)
  {
    seedToAddTo->insert_cluster_key(*citer);
  }
}
//____________________________________________________________________________..
int PHCosmicSeedCombiner::End(PHCompositeNode*)
{
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
