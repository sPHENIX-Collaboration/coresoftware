
#include "PHCosmicTrackMerger.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrClusterContainer.h>

#include <trackbase_historic/TrackSeedContainer.h>
namespace
{
  template <class T>
  inline constexpr T square(T &x)
  {
    return x * x;
  }
  template <class T>
  inline constexpr T r(T &x, T &y)
  {
    return std::sqrt(square(x) + square(y));
  }

}  // namespace
//____________________________________________________________________________..
PHCosmicTrackMerger::PHCosmicTrackMerger(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
PHCosmicTrackMerger::~PHCosmicTrackMerger()
{
}

//____________________________________________________________________________..
int PHCosmicTrackMerger::Init(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHCosmicTrackMerger::InitRun(PHCompositeNode *topNode)
{
  m_seeds = findNode::getClass<TrackSeedContainer>(topNode, "SvtxTrackSeedContainer");
  if (!m_seeds)
  {
    std::cout << PHWHERE << "No track seed container, cannot continue" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_siliconSeeds = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  m_tpcSeeds = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!m_siliconSeeds or !m_tpcSeeds)
  {
    std::cout << PHWHERE << "Missing seed container, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_geometry)
  {
    std::cout << PHWHERE << "no acts geometry, can't continue" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterMap)
  {
    std::cout << PHWHERE << "no cluster map, can't continue" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHCosmicTrackMerger::process_event(PHCompositeNode *)
{
  for (auto tr1it = m_seeds->begin(); tr1it != m_seeds->end();
       ++tr1it)
  {
    auto track1 = *tr1it;
    if (!track1)
    {
      continue;
    }
    unsigned int tpcid1 = track1->get_tpc_seed_index();
    unsigned int siid1 = track1->get_silicon_seed_index();
    auto tpcseed1 = m_tpcSeeds->get(tpcid1);
    auto silseed1 = m_siliconSeeds->get(siid1);

    TrackFitUtils::position_vector_t tr1_rz_pts, tr1_xy_pts;
    auto globTr1 = getGlobalPositions(tpcseed1);
    for (auto &pos : globTr1.second)
    {
      float clusr = r(pos.x(), pos.y());
      if (pos.y() < 0) clusr *= -1;
      tr1_rz_pts.push_back(std::make_pair(pos.z(), clusr));
      tr1_xy_pts.push_back(std::make_pair(pos.x(), pos.y()));
    }

    auto xyTr1Params = TrackFitUtils::line_fit(tr1_xy_pts);
    auto rzTr1Params = TrackFitUtils::line_fit(tr1_rz_pts);
    // float tr1xyint = std::get<1>(xyTr1Params);
    float tr1xyslope = std::get<0>(xyTr1Params);
    // float tr1rzint = std::get<1>(rzTr1Params);
    float tr1rzslope = std::get<0>(rzTr1Params);
    //! Check if the rz slope is close to 0 corresponding to an chain of clusters
    //! from an ion tail
    if (fabs(tr1rzslope) < 0.005)
    {
      m_seeds->erase(m_seeds->index(tr1it));
      continue;
    }

    for (auto tr2it = tr1it; tr2it != m_seeds->end();
         ++tr2it)
    {
      if (tr1it == tr2it)
      {
        continue;
      }

      auto track2 = *tr2it;
      if (!track2)
      {
        continue;
      }
      unsigned int tpcid2 = track2->get_tpc_seed_index();
      unsigned int siid2 = track2->get_silicon_seed_index();
      auto tpcseed2 = m_tpcSeeds->get(tpcid2);
      auto silseed2 = m_siliconSeeds->get(siid2);

      TrackFitUtils::position_vector_t tr2_rz_pts, tr2_xy_pts;
      auto globTr2 = getGlobalPositions(tpcseed2);

      for (auto &pos : globTr2.second)
      {
        tr2_rz_pts.push_back(std::make_pair(pos.z(), r(pos.x(), pos.y())));
        tr2_xy_pts.push_back(std::make_pair(pos.x(), pos.y()));
      }

      auto xyTr2Params = TrackFitUtils::line_fit(tr2_xy_pts);
      auto rzTr2Params = TrackFitUtils::line_fit(tr2_rz_pts);
      // float tr2xyint = std::get<1>(xyTr2Params);
      float tr2xyslope = std::get<0>(xyTr2Params);
      // float tr2rzint = std::get<1>(rzTr2Params);
      float tr2rzslope = std::get<0>(rzTr2Params);

      std::vector<TrkrDefs::cluskey> ckeyUnion;
      std::set_intersection(globTr1.first.begin(), globTr1.first.end(),
                            globTr2.first.begin(), globTr2.first.end(), std::back_inserter(ckeyUnion));
      if (
          //! check on common cluskeys
          (ckeyUnion.size() > 10) or
          //! check if xy/rz line fits are similar
          (fabs(tr1xyslope - tr2xyslope) < 0.5 &&
           //! rz line fits are swapped in sign because they are WRT (0,0,0)
           fabs(tr1rzslope - tr2rzslope * -1) < 0.5))
      {
        if (Verbosity() > 3)
        {
          std::cout << "Combining tr" << m_seeds->index(tr1it) << " with tr "
                    << m_seeds->index(tr2it) << " with slopes "
                    << tr1xyslope << ", " << tr2xyslope << ", "
                    << tr1rzslope << ", " << tr2rzslope << " and ckey union "
                    << ckeyUnion.size() << std::endl;
        }
        addKeys(tpcseed1, tpcseed2);
        addKeys(silseed1, silseed2);
        //! erase track 2
        m_seeds->erase(m_seeds->index(tr2it));
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHCosmicTrackMerger::addKeys(TrackSeed *toAddTo, TrackSeed *toAdd)
{
  for (auto it = toAdd->begin_cluster_keys(); it != toAdd->end_cluster_keys();
       ++it)
  {
    if (Verbosity() > 3)
    {
      auto clus = m_clusterMap->findCluster(*it);
      auto glob = m_geometry->getGlobalPosition(*it, clus);
      std::cout << "adding " << *it << " with pos " << glob.transpose() << std::endl;
    }
    toAddTo->insert_cluster_key(*it);
  }
}
PHCosmicTrackMerger::KeyPosMap
PHCosmicTrackMerger::getGlobalPositions(TrackSeed *seed)
{
  KeyPosMap glob;
  std::vector<TrkrDefs::cluskey> ckeys;
  std::vector<Acts::Vector3> positions;
  for (auto it = seed->begin_cluster_keys(); it != seed->end_cluster_keys();
       ++it)
  {
    auto key = *it;
    auto clus = m_clusterMap->findCluster(key);
    ckeys.push_back(key);
    positions.push_back(m_geometry->getGlobalPosition(key, clus));
  }
  glob.first = ckeys;
  glob.second = positions;
  return glob;
}
//____________________________________________________________________________..
int PHCosmicTrackMerger::End(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
