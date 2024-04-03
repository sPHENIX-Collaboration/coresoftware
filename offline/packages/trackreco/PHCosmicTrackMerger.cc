
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
  if (Verbosity() > 3)
  {
    std::cout << "Seed container size " << m_seeds->size() << std::endl;
  }
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

    for (auto tr2it = tr1it; tr2it != m_seeds->end();
         ++tr2it)
    {
      //! update track 1 in case more clusters have been added
      TrackFitUtils::position_vector_t tr1_rz_pts, tr1_xy_pts;
      auto globTr1 = getGlobalPositions(tpcseed1);
      for (auto &pos : globTr1.second)
      {
        float clusr = r(pos.x(), pos.y());
        if (pos.y() < 0) clusr *= -1;
        tr1_rz_pts.push_back(std::make_pair(pos.z(), clusr));
        tr1_xy_pts.push_back(std::make_pair(pos.x(), pos.y()));
      }

      float tr1xyslope = NAN;
      if (m_zeroField)
      {
        auto xyTr1Params = TrackFitUtils::line_fit(tr1_xy_pts);
        tr1xyslope = std::get<0>(xyTr1Params);
      }
      else
      {
        //! If field on, treat the circle radius as "slope"
        auto xyTr1Params = TrackFitUtils::circle_fit_by_taubin(tr1_xy_pts);
        tr1xyslope = std::get<0>(xyTr1Params);
      }

      auto rzTr1Params = TrackFitUtils::line_fit(tr1_rz_pts);
      float tr1rzslope = std::get<0>(rzTr1Params);
      //! Check if the rz slope is close to 0 corresponding to an chain of clusters
      //! from an ion tail
      if (fabs(tr1rzslope) < 0.005)
      {
        m_seeds->erase(m_seeds->index(tr1it));
        break;
      }
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
        float clusr = r(pos.x(), pos.y());
        if (pos.y() < 0) clusr *= -1;
        tr2_rz_pts.push_back(std::make_pair(pos.z(), clusr));
        tr2_xy_pts.push_back(std::make_pair(pos.x(), pos.y()));
      }

      float tr2xyslope = NAN;
      if (m_zeroField)
      {
        auto xyTr2Params = TrackFitUtils::line_fit(tr2_xy_pts);
        tr2xyslope = std::get<0>(xyTr2Params);
      }
      else
      {
        //! If field on, treat the circle radius as "slope"
        auto xyTr2Params = TrackFitUtils::circle_fit_by_taubin(tr2_xy_pts);
        tr2xyslope = std::get<0>(xyTr2Params);
      }

      auto rzTr2Params = TrackFitUtils::line_fit(tr2_rz_pts);
      // float tr2xyint = std::get<1>(xyTr2Params);
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
          std::cout << "Combining tr" << m_seeds->index(tr1it) << " with sil/tpc " << tpcid1
                    << ", " << siid1 << " with tr "
                    << m_seeds->index(tr2it) << " with sil/tpc " << tpcid2 << ", "
                    << siid2 << " with slopes "
                    << tr1xyslope << ", " << tr2xyslope << ", "
                    << tr1rzslope << ", " << tr2rzslope << " and ckey union "
                    << ckeyUnion.size() << std::endl;
        }

        for (auto &key : globTr2.first)
        {
          globTr1.first.push_back(key);
        }
        for (auto &pos : globTr2.second)
        {
          globTr1.second.push_back(pos);
        }

        addKeys(tpcseed1, tpcseed2);
        if (silseed1 && silseed2)
        {
          addKeys(silseed1, silseed2);
        }
        //! erase track 2
        m_seeds->erase(m_seeds->index(tr2it));
      }
    }

    //! remove any obvious outlier clusters from the track that were mistakenly
    //! picked up by the seeder
    if (m_iter == 1)
    {
      getBestClustersPerLayer(tpcseed1);
      if (silseed1)
      {
        getBestClustersPerLayer(silseed1);
      }
    }
    removeOutliers(tpcseed1);
    if (silseed1)
    {
      removeOutliers(silseed1);
    }
  }
  if (Verbosity() > 3)
  {
    std::cout << "Seed container size to finish " << m_seeds->size() << std::endl;
    for (auto &seed : *m_seeds)
    {
      if (seed)
      {
        seed->identify();
        std::cout << std::endl;
        unsigned int tpcid1 = seed->get_tpc_seed_index();
        unsigned int siid1 = seed->get_silicon_seed_index();
        auto tpcseed1 = m_tpcSeeds->get(tpcid1);
        auto silseed1 = m_siliconSeeds->get(siid1);
        tpcseed1->identify();
        if (silseed1)
        {
          silseed1->identify();
        }
      }
      else
      {
        std::cout << "nullptr seed was removed" << std::endl;
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHCosmicTrackMerger::getBestClustersPerLayer(TrackSeed *seed)
{
  TrackFitUtils::position_vector_t tr_rz_pts, tr_xy_pts;

  auto glob = getGlobalPositions(seed);

  for (const auto &pos : glob.second)
  {
    float clusr = r(pos.x(), pos.y());
    if (pos.y() < 0)
    {
      clusr *= -1;
    }
    // skip tpot clusters, as they are always bad in 1D due to 1D resolution
    if (fabs(clusr) > 80.)
    {
      continue;
    }
    tr_rz_pts.push_back(std::make_pair(pos.z(), clusr));
    tr_xy_pts.push_back(std::make_pair(pos.x(), pos.y()));
  }

  auto xyParams = TrackFitUtils::line_fit(tr_xy_pts);
  auto rzParams = TrackFitUtils::line_fit(tr_rz_pts);
  std::map<int, std::pair<float, float>> bestLayerDcasxy, bestLayerDcasrz;
  std::map<int, std::pair<TrkrDefs::cluskey, TrkrDefs::cluskey>> bestLayerCluskeys;
  for (int i = 0; i < 57; i++)
  {
    bestLayerDcasrz.insert(std::make_pair(i, std::make_pair(std::numeric_limits<float>::max(), std::numeric_limits<float>::max())));
    bestLayerDcasxy.insert(std::make_pair(i, std::make_pair(std::numeric_limits<float>::max(), std::numeric_limits<float>::max())));
    bestLayerCluskeys.insert(std::make_pair(i, std::make_pair(0, 0)));
  }

  std::vector<TrkrDefs::cluskey> tpotClus;
  for (int i = 0; i < glob.first.size(); i++)
  {
    auto &pos = glob.second[i];
    float clusr = r(pos.x(), pos.y());
    if (pos.y() < 0)
    {
      clusr *= -1;
    }

    float perpxyslope = -1. / std::get<0>(xyParams);
    float perpxyint = pos.y() - perpxyslope * pos.x();
    float perprzslope = -1. / std::get<0>(rzParams);
    float perprzint = clusr - perprzslope * pos.z();

    float pcax = (perpxyint - std::get<1>(xyParams)) / (std::get<0>(xyParams) - perpxyslope);
    float pcay = std::get<0>(xyParams) * pcax + std::get<1>(xyParams);

    float pcaz = (perprzint - std::get<1>(rzParams)) / (std::get<0>(rzParams) - perprzslope);
    float pcar = std::get<0>(rzParams) * pcaz + std::get<1>(rzParams);
    float dcax = pcax - pos.x();
    float dcay = pcay - pos.y();
    float dcar = pcar - clusr;
    float dcaz = pcaz - pos.z();
    float dcaxy = std::sqrt(square(dcax) + square(dcay));
    float dcarz = std::sqrt(square(dcar) + square(dcaz));
    auto trkid = TrkrDefs::getTrkrId(glob.first[i]);
    auto layer = TrkrDefs::getLayer(glob.first[i]);
    //! combine layers 3/4 and 5/6 since they are half layers
    if (layer == 4)
    {
      layer = 3;
    }
    if (layer == 6)
    {
      layer = 5;
    }
    float dcaxydiff1 = std::fabs(dcaxy - bestLayerDcasxy.find(layer)->second.first);
    float dcaxydiff2 = std::fabs(dcaxy - bestLayerDcasxy.find(layer)->second.second);
    if (trkid == TrkrDefs::TrkrId::mvtxId || trkid == TrkrDefs::TrkrId::tpcId)
    {
      if (dcaxydiff1 > dcaxydiff2)
      {
        if (dcaxy < bestLayerDcasxy.find(layer)->second.first &&
            dcarz < bestLayerDcasrz.find(layer)->second.first)
        {
          bestLayerDcasxy.find(layer)->second.first = dcaxy;
          bestLayerDcasrz.find(layer)->second.first = dcarz;
          bestLayerCluskeys.find(layer)->second.first = glob.first[i];
        }
      }
      else
      {
        if (dcaxy < bestLayerDcasxy.find(layer)->second.second &&
            dcarz < bestLayerDcasrz.find(layer)->second.second)
        {
          bestLayerDcasxy.find(layer)->second.second = dcaxy;
          bestLayerDcasrz.find(layer)->second.second = dcarz;
          bestLayerCluskeys.find(layer)->second.second = glob.first[i];
        }
      }
    }
    else if (trkid == TrkrDefs::TrkrId::inttId)
    {
      if (dcaxydiff1 > dcaxydiff2)
      {
        if (dcaxy < bestLayerDcasxy.find(layer)->second.first)
        {
          bestLayerDcasxy.find(layer)->second.first = dcaxy;
          bestLayerDcasrz.find(layer)->second.first = dcarz;
          bestLayerCluskeys.find(layer)->second.first = glob.first[i];
        }
      }
      else
      {
        if (dcaxy < bestLayerDcasxy.find(layer)->second.second)
        {
          bestLayerDcasxy.find(layer)->second.second = dcaxy;
          bestLayerDcasrz.find(layer)->second.second = dcarz;
          bestLayerCluskeys.find(layer)->second.second = glob.first[i];
        }
      }
    }
    else if (trkid == TrkrDefs::TrkrId::micromegasId)
    {
      //! add all tpot clusters, since they are 1d
      tpotClus.emplace_back(glob.first[i]);
    }
  }
  seed->identify();
  // now erase all cluskeys and fill with the new set of cluskeys
  seed->clear_cluster_keys();
  for (auto &[layer, keypair] : bestLayerCluskeys)
  {
    for (auto &key : {keypair.first, keypair.second})
      if (key > 0)
      {
        seed->insert_cluster_key(key);
      }
  }
  for (auto &key : tpotClus)
  {
    seed->insert_cluster_key(key);
  }
}
void PHCosmicTrackMerger::removeOutliers(TrackSeed *seed)
{
  TrackFitUtils::position_vector_t tr_rz_pts, tr_xy_pts;
  auto glob = getGlobalPositions(seed);
  for (const auto &pos : glob.second)
  {
    float clusr = r(pos.x(), pos.y());
    if (pos.y() < 0) clusr *= -1;
    // skip tpot clusters, as they are always bad in 1D due to 1D resolution
    if (fabs(clusr) > 80.) continue;
    tr_rz_pts.push_back(std::make_pair(pos.z(), clusr));
    tr_xy_pts.push_back(std::make_pair(pos.x(), pos.y()));
  }

  auto xyParams = TrackFitUtils::line_fit(tr_xy_pts);
  auto rzParams = TrackFitUtils::line_fit(tr_rz_pts);

  for (int i = 0; i < glob.first.size(); i++)
  {
    auto &pos = glob.second[i];
    float clusr = r(pos.x(), pos.y());
    if (pos.y() < 0)
    {
      clusr *= -1;
    }
    // skip tpot clusters, as they are always bad in 1D due to 1D resolution
    if (fabs(clusr) > 80.)
    {
      continue;
    }
    float perpxyslope = -1. / std::get<0>(xyParams);
    float perpxyint = pos.y() - perpxyslope * pos.x();
    float perprzslope = -1. / std::get<0>(rzParams);
    float perprzint = clusr - perprzslope * pos.z();

    float pcax = (perpxyint - std::get<1>(xyParams)) / (std::get<0>(xyParams) - perpxyslope);
    float pcay = std::get<0>(xyParams) * pcax + std::get<1>(xyParams);

    float pcaz = (perprzint - std::get<1>(rzParams)) / (std::get<0>(rzParams) - perprzslope);
    float pcar = std::get<0>(rzParams) * pcaz + std::get<1>(rzParams);
    float dcax = pcax - pos.x();
    float dcay = pcay - pos.y();
    float dcar = pcar - clusr;
    float dcaz = pcaz - pos.z();
    float dcaxy = std::sqrt(square(dcax) + square(dcay));
    float dcarz = std::sqrt(square(dcar) + square(dcaz));

    //! exclude silicon from DCAz cut
    if (dcaxy > m_dcaxycut ||
        (dcarz > m_dcarzcut && (r(pos.x(), pos.y()) > 20)))
    {
      if (Verbosity() > 2)
      {
        std::cout << "Erasing ckey " << glob.first[i] << " with position " << pos.transpose() << std::endl;
      }
      seed->erase_cluster_key(glob.first[i]);
    }
  }
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
