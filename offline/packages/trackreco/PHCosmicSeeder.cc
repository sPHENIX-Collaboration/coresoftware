
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

#include <TFile.h>
#include <TNtuple.h>
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

  if(m_analysis)
  {
  m_outfile = new TFile("PHCosmicSeeder.root", "recreate");
  m_tup = new TNtuple("ntp_seed", "seed", "event:seed:nclus:xyint:xyslope:rzint:rzslope");
  }

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
  if (Verbosity() > 2)
  {
    std::cout << "cluster map size is " << clusterPositions.size() << std::endl;
  }
  auto seeds = makeSeeds(clusterPositions);

  std::sort(seeds.begin(), seeds.end(),
            [](seed a, seed b)
            { return a.ckeys.size() > b.ckeys.size(); });

  auto prunedSeeds = combineSeeds(seeds, clusterPositions);
  int iseed = 0;
  std::cout << "pruned seed size " << prunedSeeds.size() << std::endl;

  for (auto& seed : prunedSeeds)
  {
 
    std::cout << "seed qualities "
              << seed.xyslope << ", " << seed.xyintercept << ", " << seed.rzslope << ", "
              << seed.rzintercept << std::endl;
              if(m_analysis)
              {
                float seed_data[] = {
                  (float) m_event,
                  (float) iseed,
                  (float) seed.ckeys.size(),
                  seed.xyintercept,
                  seed.xyslope,
                  seed.rzintercept,
                  seed.rzslope
                };
                m_tup->Fill(seed_data);
              }
    auto svtxseed = std::make_unique<TrackSeed_v1>();
    for (auto& key : seed.ckeys)
    {
      svtxseed->insert_cluster_key(key);
    }
    std::cout << "seed size " << svtxseed->size_cluster_keys() << std::endl;
    m_seedContainer->insert(svtxseed.get());
    ++iseed;
  }
  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}
PHCosmicSeeder::SeedVector PHCosmicSeeder::combineSeeds(PHCosmicSeeder::SeedVector& initialSeeds,
                                                        PHCosmicSeeder::PositionMap& clusterPositions)
{
  SeedVector prunedSeeds;
  std::set<int> seedsToDelete;

  for (int i = 0; i < initialSeeds.size(); ++i)
  {
    for (int j = i; j < initialSeeds.size(); ++j)
    {
      if (i == j)
      {
        continue;
      }
      auto seed1 = initialSeeds[i];
      auto seed2 = initialSeeds[j];

      // recalculate seed parameters
      recalculateSeedLineParameters(seed1, clusterPositions, true);
      recalculateSeedLineParameters(seed2, clusterPositions, true);
      recalculateSeedLineParameters(seed1, clusterPositions, false);
      recalculateSeedLineParameters(seed2, clusterPositions, false);
      std::cout << "seed after qualities "
                << seed1.xyslope << ", " << seed1.xyintercept << ", " << seed1.rzslope << ", "
                << seed1.rzintercept << std::endl;
      std::cout << "seed after qualities "
                << seed2.xyslope << ", " << seed2.xyintercept << ", " << seed2.rzslope << ", "
                << seed2.rzintercept << std::endl;
      if (fabs(seed1.xyslope - seed2.xyslope) < 0.1 &&
          fabs(seed1.xyintercept - seed2.xyintercept) < 2. &&
          fabs(seed1.rzslope - seed2.rzslope) < 0.1 &&
          fabs(seed1.rzintercept - seed2.rzintercept) < 2.)
      {
        for (auto& key : seed2.ckeys)
        {
          seed1.ckeys.insert(key);
        }
        seedsToDelete.insert(j);
      }
    }
  }

  for (int i = 0; i < initialSeeds.size(); ++i)
  {
    if (seedsToDelete.find(i) != seedsToDelete.end() or initialSeeds[i].ckeys.size() < 5)
    {
      continue;
    }
    prunedSeeds.push_back(initialSeeds[i]);
  }

  return prunedSeeds;
}
PHCosmicSeeder::SeedVector
PHCosmicSeeder::makeSeeds(PHCosmicSeeder::PositionMap& clusterPositions)
{
  PHCosmicSeeder::SeedVector seeds;
  std::set<TrkrDefs::cluskey> keys;
  int seednum = 0;
  for (auto& [key1, pos1] : clusterPositions)
  {
    for (auto& [key2, pos2] : clusterPositions)
    {
      if (key1 == key2)
      {
        continue;
      }
      // make a cut on clusters to at least be close to each other within a few cm
      float dist = (pos2 - pos1).norm();
      if (Verbosity() > 5)
      {
        std::cout << "checking keys " << key1 << ", " << key2 << " with dist apart " << dist << std::endl;
        std::cout << "positions are " << pos1.transpose() << "    ,     "
                  << pos2.transpose() << std::endl;
      }
      if (dist > 2.)
      {
        continue;
      }
      PHCosmicSeeder::seed doub;
      // std::cout << "pos 1 and 2 " << pos1.transpose() << " , " << pos2.transpose() << std::endl;
      doub.xyslope = (pos2.y() - pos1.y()) / (pos2.x() - pos1.x());
      doub.xyintercept = pos1.y() - doub.xyslope * pos1.x();
      doub.rzslope = (r(pos2.x(), pos2.y()) - r(pos1.x(), pos1.y())) / (pos2.z() - pos1.z());
      doub.rzintercept = pos1.z() * doub.rzslope + r(pos1.x(), pos1.y());
      // std::cout << "vals are " << doub.rzslope << ", " << doub.rzintercept << std::endl;
      keys.insert(key1);
      keys.insert(key2);
      doub.ckeys = keys;
      keys.clear();
      seeds.push_back(doub);
      seednum++;
    }
  }
  if (Verbosity() > 2)
  {
    std::cout << "odublet sizes " << seeds.size() << std::endl;
  }
  for (auto& dub : seeds)
  {
    if (Verbosity() > 2)
    {
      std::cout << "doublet has " << dub.ckeys.size() << " keys " << std::endl;
      for (auto key : dub.ckeys)
      {
        std::cout << "position is " << clusterPositions.find(key)->second.transpose() << std::endl;
      }
    }
    auto begin = dub.ckeys.begin();
    auto pos1 = clusterPositions.find(*(begin))->second;
    std::advance(begin, 1);
    auto pos2 = clusterPositions.find(*(begin))->second;

    for (auto& [key, pos] : clusterPositions)
    {
      //! skip existing keys
      if (std::find(dub.ckeys.begin(), dub.ckeys.end(), key) != dub.ckeys.end())
      {
        continue;
      }
      // only look at the cluster that is within 2cm of the doublet clusters
      float dist1 = (pos1 - pos).norm();
      float dist2 = (pos2 - pos).norm();
      if (dist1 < dist2)
      {
        if (dist1 > 2.)
          continue;
      }
      else
      {
        if (dist2 > 2.)
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
        if (Verbosity() > 2)
        {
          std::cout << "   adding ckey " << key << " with box dca "
                    << predy << ", " << pos.y() << " and " << predr << ", " << r(pos.x(), pos.y())
                    << std::endl;
        }
        dub.ckeys.insert(key);
      }
    }
  }

  /// erase doublets
  for (auto it = seeds.begin(); it != seeds.end(); ++it)
  {
    auto seed = *it;
    std::cout << "double seed size " << seed.ckeys.size() << std::endl;
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
void PHCosmicSeeder::recalculateSeedLineParameters(seed& seed,
                                               PHCosmicSeeder::PositionMap& clusters, bool isXY)
{
  float avgx = 0;
  float avgy = 0;
  PHCosmicSeeder::PositionMap seedClusters;
  for (auto& key : seed.ckeys)
  {
    auto glob = clusters.find(key)->second;
    if (isXY)
    {
      avgx += glob.x();
      avgy += glob.y();
    }
    else
    {
      avgx += glob.z();
      avgy += r(glob.x(), glob.y());
    }
    seedClusters.insert(std::make_pair(key, glob));
  }

  avgx /= seed.ckeys.size();
  avgy /= seed.ckeys.size();
  float num = 0;
  float denom = 0;
  for (auto& [key, glob] : seedClusters)
  {
    if (isXY)
    {
      num += (glob.x() - avgx) * (glob.y() - avgy);
      denom += square(glob.x() - avgx);
    }
    else
    {
      num += (glob.z() - avgx) * (r(glob.x(), glob.y()) - avgy);
      denom += square(glob.z() - avgx);
    }
  }
  if (isXY)
  {
    seed.xyslope = num / denom;
    seed.xyintercept = avgy - seed.xyslope * avgx;
  }
  else
  {
    seed.rzslope = num / denom;
    seed.rzintercept = avgy - seed.rzslope * avgx;
  }
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
  if(m_outfile)
  {
  m_outfile->cd();
  m_tup->Write();
  m_outfile->Close();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
