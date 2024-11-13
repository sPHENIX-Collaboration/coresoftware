
#include "PHSiliconCosmicSeeding.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v2.h>
#include <trackbase_historic/TrackSeedHelper.h>

namespace
{
  template <class T>
  inline constexpr T square(const T &x)
  {
    return x * x;
  }
  template <class T>
  inline constexpr T r(const T &x, const T &y)
  {
    return std::sqrt(square(x) + square(y));
  }
}  // namespace

//____________________________________________________________________________..
PHSiliconCosmicSeeding::PHSiliconCosmicSeeding(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
PHSiliconCosmicSeeding::~PHSiliconCosmicSeeding() = default;

//____________________________________________________________________________..
int PHSiliconCosmicSeeding::Init(PHCompositeNode * /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHSiliconCosmicSeeding::InitRun(PHCompositeNode *topNode)
{
  int ret = createNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }
  ret = getNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHSiliconCosmicSeeding::process_event(PHCompositeNode * /*unused*/)
{
  PHSiliconCosmicSeeding::PositionMap clusterPositions;
  std::set<TrkrDefs::TrkrId> detectors;
  detectors.insert(TrkrDefs::TrkrId::mvtxId);
  detectors.insert(TrkrDefs::TrkrId::inttId);
  for (const auto &det : detectors)
  {
    for (const auto &hitsetkey : m_clusterContainer->getHitSetKeys(det))
    {
      auto range = m_clusterContainer->getClusters(hitsetkey);
      for (auto citer = range.first; citer != range.second; ++citer)
      {
        const auto ckey = citer->first;
        const auto cluster = citer->second;
        const auto global = m_tGeometry->getGlobalPosition(ckey, cluster);
        clusterPositions.insert(std::make_pair(ckey, global));
      }
    }
  }
  if (Verbosity() > 2)
  {
    std::cout << "cluster map size is " << clusterPositions.size() << std::endl;
  }

  //! to protect against events with hot channels. Cosmics should not produce more than
  //! 500 clusters in the silicon
  if (clusterPositions.size() > 500)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  auto doublets = makeDoublets(clusterPositions);

  if (Verbosity() > 2)
  {
    std::cout << "doublets size is " << doublets.size() << std::endl;
    if (doublets.size() > 0)
    {
      std::cout << "nonzero doublet size" << std::endl;
    }
  }

  auto longseeds = addClustersOnLine(doublets, clusterPositions);

  if (Verbosity() > 2)
  {
    std::cout << "longseeds size is " << longseeds.size() << std::endl;
  }
  auto finalseeds = combineSeeds(longseeds);
  if (Verbosity() > 2)
  {
    std::cout << "final seeds size is " << finalseeds.size() << std::endl;
  }

  pruneSeeds(finalseeds, clusterPositions);
  for (auto &s : finalseeds)
  {
    //! make some quality cuts on the seeds
    if (s.ckeys.size() < 3)
    {
      continue;
    }
    int nmaps = 0;
    int nintt = 0;
    for (auto &key : s.ckeys)
    {
      if (TrkrDefs::getTrkrId(key) == TrkrDefs::TrkrId::mvtxId)
      {
        nmaps++;
      }
      if (TrkrDefs::getTrkrId(key) == TrkrDefs::TrkrId::inttId)
      {
        nintt++;
      }
    }
    if (nmaps > 3 && nmaps < 9 && nintt > 2 && nintt < 7)
    {
      std::unique_ptr<TrackSeed_v2> si_seed = std::make_unique<TrackSeed_v2>();
      for (auto &key : s.ckeys)
      {
        si_seed->insert_cluster_key(key);
      }
      TrackSeedHelper::circleFitByTaubin(si_seed.get(), clusterPositions, 0, 7);
      TrackSeedHelper::lineFit(si_seed.get(), clusterPositions, 0, 7);

      // calculate phi and assign
      si_seed->set_phi(TrackSeedHelper::get_phi(si_seed.get(), clusterPositions));

      if (Verbosity() > 3)
      {
        si_seed->identify();
      }
      m_seedContainer->insert(si_seed.get());
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHSiliconCosmicSeeding::pruneSeeds(SeedVector &seeds, PositionMap &clusterPositions)
{
  SeedVector prunedSeeds;
  for (auto &s : seeds)
  {
    TrackFitUtils::position_vector_t xypoints;
    for (auto &key : s.ckeys)
    {
      auto pos = clusterPositions.find(key)->second;
      xypoints.push_back(std::make_pair(pos.x(), pos.y()));
    }

    auto xyLineParams = TrackFitUtils::line_fit(xypoints);
    float lineSlope = std::get<0>(xyLineParams);
    float lineIntercept = std::get<1>(xyLineParams);
    std::set<TrkrDefs::cluskey> newKeys;
    for (auto &key : s.ckeys)
    {
      auto pos = clusterPositions.find(key)->second;
      float distance = std::abs(lineSlope * pos.x() - pos.y() + lineIntercept) / std::sqrt(lineSlope * lineSlope + 1);
      if (distance < 0.5)
      {
        newKeys.insert(key);
      }
    }
    s.ckeys = newKeys;
  }
  return;
}
PHSiliconCosmicSeeding::SeedVector
PHSiliconCosmicSeeding::combineSeeds(SeedVector &seeds)
{
  SeedVector newseeds;
  std::vector<TrkrDefs::cluskey> overlapKeys;
  for (size_t i = 0; i < seeds.size(); i++)
  {
    auto &s1 = seeds[i];
    for (size_t j = i; j < seeds.size(); j++)
    {
      auto &s2 = seeds[j];

      std::set_intersection(s1.ckeys.begin(), s1.ckeys.end(),
                            s2.ckeys.begin(), s2.ckeys.end(), std::back_inserter(overlapKeys));
      if (overlapKeys.size() > 1)
      {
        if (Verbosity() > 3)
        {
          std::cout << "duplicate seeds found with " << std::endl;
          for (auto &key : s1.ckeys)
          {
            std::cout << key << ", ";
          }
          std::cout << std::endl;
          for (auto &key : s2.ckeys)
          {
            std::cout << key << ", ";
          }
          std::cout << std::endl;
        }
        //! check it isn't already in there
        bool found = false;
        for (auto &newseed : newseeds)
        {
          overlapKeys.clear();
          std::set_intersection(s1.ckeys.begin(), s1.ckeys.end(),
                                newseed.ckeys.begin(), newseed.ckeys.end(), std::back_inserter(overlapKeys));
          if (overlapKeys.size() > 1)
          {
            found = true;
            break;
          }
        }
        if (found)
        {
          break;
        }

        newseeds.push_back(s1);
        break;
      }
      overlapKeys.clear();
    }
  }
  if (Verbosity() > 3)
  {
    std::cout << "final seeds " << std::endl;
    for (auto &s : newseeds)
    {
      for (auto &key : s.ckeys)
      {
        std::cout << key << ", ";
      }
      std::cout << std::endl;
    }
  }
  return newseeds;
}
PHSiliconCosmicSeeding::SeedVector PHSiliconCosmicSeeding::addClustersOnLine(SeedVector &doublets, PositionMap &clusterPositions)
{
  SeedVector longseeds;
  for (auto doublet : doublets)
  {
    TrackFitUtils::position_vector_t xypoints;
    for (auto &key : doublet.ckeys)
    {
      auto pos = clusterPositions.find(key)->second;
      xypoints.push_back(std::make_pair(pos.x(), pos.y()));
    }
    std::vector<TrkrDefs::cluskey> newClusKeysxy;
    std::vector<Acts::Vector3> newClusPosxy;
    auto xyparams = TrackFitUtils::line_fit(xypoints);
    float nxyclusters = TrackFitUtils::addClustersOnLine(xyparams, true, 0.7,
                                                         m_tGeometry, m_clusterContainer, newClusPosxy, newClusKeysxy, 0, 7);
    if (nxyclusters > 1)
    {
      for (auto &newkey : newClusKeysxy)
      {
        doublet.ckeys.insert(newkey);
      }
    }
    if (doublet.ckeys.size() > 3 && doublet.ckeys.size() < 20)
    {
      longseeds.push_back(doublet);
    }
  }

  if (Verbosity() > 2)
  {
    std::cout << "num seeds " << longseeds.size() << std::endl;
    int i = 0;
    for (auto &s : longseeds)
    {
      std::cout << "seed has " << s.ckeys.size() << " clusters" << std::endl;
      std::cout << "seed " << i << " has keys " << std::endl;
      for (auto &key : s.ckeys)
      {
        std::cout << "    " << key << ", ";
      }
      std::cout << std::endl;
      i++;
    }
  }

  return longseeds;
}
PHSiliconCosmicSeeding::SeedVector PHSiliconCosmicSeeding::makeDoublets(PositionMap &clusterPositions)
{
  SeedVector doublets;
  std::set<TrkrDefs::cluskey> keys;

  for (auto iter = clusterPositions.begin(); iter != clusterPositions.end(); ++iter)
  {
    auto key1 = iter->first;
    auto pos1 = iter->second;
    for (auto iter2 = iter; iter2 != clusterPositions.end(); ++iter2)
    {
      auto key2 = iter2->first;
      auto pos2 = iter2->second;
      if (key1 == key2)
      {
        continue;
      }
      float dist = (pos2 - pos1).norm();
      if (Verbosity() > 5)
      {
        std::cout << "checking keys " << key1 << ", " << key2 << " with dist apart " << dist << std::endl;
        std::cout << "positions are " << pos1.transpose() << "    ,     "
                  << pos2.transpose() << std::endl;
      }
      if (dist > m_maxDoubletDistance)
      {
        continue;
      }
      // skip doublets that are too close together
      if (dist < 0.1)
      {
        continue;
      }
      PHSiliconCosmicSeeding::seed doub;
      doub.xyslope = (pos2.y() - pos1.y()) / (pos2.x() - pos1.x());
      doub.xyintercept = pos1.y() - doub.xyslope * pos1.x();
      doub.rzslope = (r(pos2.x(), pos2.y()) - r(pos1.x(), pos1.y())) / (pos2.z() - pos1.z());
      doub.rzintercept = pos1.z() * doub.rzslope + r(pos1.x(), pos1.y());

      keys.insert(key1);
      keys.insert(key2);
      doub.ckeys = keys;
      keys.clear();
      doublets.push_back(doub);
    }
  }
  return doublets;
}
//____________________________________________________________________________..
int PHSiliconCosmicSeeding::End(PHCompositeNode * /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSiliconCosmicSeeding::getNodes(PHCompositeNode *topNode)
{
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE << "No acts reco geometry, bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterContainer)
  {
    std::cout << PHWHERE << "No cluster container on node tree, bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSiliconCosmicSeeding::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHSiliconCosmicSeeding::createNodes");
  }

  PHNodeIterator dstIter(dstNode);
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_seedContainer = findNode::getClass<TrackSeedContainer>(topNode, m_trackMapName);
  if (!m_seedContainer)
  {
    m_seedContainer = new TrackSeedContainer_v1;
    PHIODataNode<PHObject> *trackNode =
        new PHIODataNode<PHObject>(m_seedContainer, m_trackMapName, "PHObject");
    svtxNode->addNode(trackNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
