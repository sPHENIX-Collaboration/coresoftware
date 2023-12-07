
#include "PHCosmicSeeder.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v1.h>

namespace
{
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }
  template <class T>
  inline constexpr T r(const T& x, const T& y)
  {
    return std::sqrt(square(x) + square(y));
  }
}  // namespace
//____________________________________________________________________________..
PHCosmicSeeder::PHCosmicSeeder(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
PHCosmicSeeder::~PHCosmicSeeder()
{
}

//____________________________________________________________________________..
int PHCosmicSeeder::Init(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHCosmicSeeder::InitRun(PHCompositeNode* topNode)
{
  int ret = getNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }
  ret = createNodes(topNode);

  return ret;
}

//____________________________________________________________________________..
int PHCosmicSeeder::process_event(PHCompositeNode*)
{
  PHCosmicSeeder::PositionMap clusterPositions;
  for (const auto& hitsetkey : m_clusterContainer->getHitSetKeys(TrkrDefs::TrkrId::tpcId))
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
  if(Verbosity() > 2)
  {
    std::cout << "cluster map size is " << clusterPositions.size() << std::endl;
  }
  auto seeds = makeSeeds(clusterPositions);

  std::sort(seeds.begin(), seeds.end(),
            [](seed a, seed b)
            { return a.ckeys.size() > b.ckeys.size(); });

  std::cout << "Total seeds found " << seeds.size() << std::endl;
  int i = 0;
  for (auto& seed : seeds)
  {
    std::cout << "processing seed " << i << std::endl;
    ++i;

    auto svtxseed = std::make_unique<TrackSeed_v1>();
    for(auto& key : seed.ckeys)
    {
      std::cout << "   " << clusterPositions.find(key)->second.transpose()
                << std::endl;
      svtxseed->insert_cluster_key(key);
    }
    m_seedContainer->insert(svtxseed.get());
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

PHCosmicSeeder::SeedVec
PHCosmicSeeder::makeSeeds(PHCosmicSeeder::PositionMap& clusterPositions)
{
  PHCosmicSeeder::SeedVec seeds;
  std::vector<TrkrDefs::cluskey> keys;
  for (auto& [key1, pos1] : clusterPositions)
  {
    for (auto& [key2, pos2] : clusterPositions)
    {
      if(key1 == key2) 
      {
        continue;
      }
      // make a cut on clusters to at least be close to each other within a few cm
      float dist = (pos2 - pos1).norm();
      if(Verbosity() > 5)
      {
        std::cout << "checking keys " << key1 << ", " << key2 << " with dist apart " << dist << std::endl;
        std::cout << "positions are " << pos1.transpose() << "    ,     "
                  << pos2.transpose() << std::endl;
      }
      if (dist > 7.)
      {
        continue;
      }
      PHCosmicSeeder::seed doub;
      //std::cout << "pos 1 and 2 " << pos1.transpose() << " , " << pos2.transpose() << std::endl;
      doub.xyslope = (pos2.y() - pos1.y()) / (pos2.x() - pos1.x());
      doub.xyintercept = pos1.y() - doub.xyslope * pos1.x();
      doub.rzslope = (r(pos2.x(), pos2.y()) - r(pos1.x(), pos1.y())) / (pos2.z() - pos1.z());
      doub.rzintercept = pos1.z() * doub.rzslope + r(pos1.x(), pos1.y());
      //std::cout << "vals are " << doub.rzslope << ", " << doub.rzintercept << std::endl;
      keys.push_back(key1);
      keys.push_back(key2);
      doub.ckeys = keys;
      keys.clear();
      seeds.push_back(doub);
    }
  }
if(Verbosity() > 2)
{
  std::cout << "odublet sizes " << seeds.size() << std::endl;
}
  for (auto& dub : seeds)
  {
    if(Verbosity() > 2)
    {
      std::cout << "doublet has " << dub.ckeys.size() << " keys " << std::endl;
      for(auto key : dub.ckeys)
      {
        std::cout << "position is " << clusterPositions.find(key)->second.transpose() << std::endl;
      }
    }
    for (auto& [key, pos] : clusterPositions)
    {
      //! skip existing keys
      if (std::find(dub.ckeys.begin(), dub.ckeys.end(), key) != dub.ckeys.end())
      {
        continue;
      }
      float predy = dub.xyslope * pos.x() + dub.xyintercept;
      float predr = dub.rzslope * pos.z() + dub.rzintercept;
      if (Verbosity() > 2)
      {
        std::cout << "testing ckey " << key << " with box dca "
                  << predy << ", " << pos.transpose() << " and " << predr << ", " << r(pos.x(), pos.y())
                  << std::endl;
      }
      if (fabs(predy - pos.y()) < m_xyTolerance) 
      {
        if(Verbosity() > 2)
        {
          std::cout << "   adding ckey " << key << " with box dca "
                    << predy << ", " << pos.y() << " and " << predr << ", " << r(pos.x(), pos.y())
                    << std::endl;
        }
        dub.ckeys.push_back(key);
      }
    }
  }

  /// erase doublets
  for (auto it = seeds.begin(); it != seeds.end(); ++it)
  {
    auto seed = *it;
    if (Verbosity() > 2)
    {
      std::cout << "keys in seed " << std::endl;
      for (auto& key : seed.ckeys)
      {
        std::cout << key << ", ";
      }
      std::cout << std::endl
                << "seed pos " << std::endl;
      for (auto& key : seed.ckeys)
      {
        std::cout << clusterPositions.find(key)->second.transpose() << std::endl;
      }
      std::cout << "Done printing seed info" << std::endl;
    }
    if (seed.ckeys.size() < 3)
    {
      seeds.erase(it);
      --it;
    }
  }

  return seeds;
}

int PHCosmicSeeder::getNodes(PHCompositeNode* topNode)
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

int PHCosmicSeeder::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHCosmicSeeder::createNodes");
  }

  PHNodeIterator dstIter(dstNode);
  PHCompositeNode* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_seedContainer = findNode::getClass<TrackSeedContainer>(topNode, m_trackMapName);
  if (!m_seedContainer)
  {
    m_seedContainer = new TrackSeedContainer_v1;
    PHIODataNode<PHObject>* trackNode =
        new PHIODataNode<PHObject>(m_seedContainer, m_trackMapName, "PHObject");
    svtxNode->addNode(trackNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
//____________________________________________________________________________..
int PHCosmicSeeder::End(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
