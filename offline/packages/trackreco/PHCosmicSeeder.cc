
#include "PHCosmicSeeder.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v2.h>
#include <trackbase/TrkrCluster.h>

#include <TFile.h>
#include <TNtuple.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <memory>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
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
    double sign = 1;
    if (y < 0)
    {
      sign = -1;
    }
    return std::sqrt(square(x) + square(y)) * sign;
  }
}  // namespace
//____________________________________________________________________________..
PHCosmicSeeder::PHCosmicSeeder(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int PHCosmicSeeder::Init(PHCompositeNode* /*unused*/)
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

  if (m_analysis)
  {
    m_outfile = new TFile("PHCosmicSeeder.root", "recreate");
    m_tup = new TNtuple("ntp_seed", "seed",
                        "event:seed:nclus:xyint:xyslope:xzint:xzslope:"
                        "longestxyint:longestxyslope:longestxzint:longestxzslope");
  }

  return ret;
}

//____________________________________________________________________________..
int PHCosmicSeeder::process_event(PHCompositeNode* /*unused*/)
{
  PHCosmicSeeder::PositionMap clusterPositions;
  for (const auto& hitsetkey : m_clusterContainer->getHitSetKeys(m_trackerId))
  {
    auto range = m_clusterContainer->getClusters(hitsetkey);
    for (auto citer = range.first; citer != range.second; ++citer)
    {
      const auto ckey = citer->first;
      const auto cluster = citer->second;
      if(cluster->getMaxAdc() < m_adcCut){
        continue;
      }
      const auto global = m_tGeometry->getGlobalPosition(ckey, cluster);
      clusterPositions.insert(std::make_pair(ckey, global));
    }
  }
  if (Verbosity() > 2)
  {
    std::cout << "cluster map size is " << clusterPositions.size() << std::endl;
  }
  auto seeds = makeSeeds(clusterPositions);

  if (Verbosity() > 1)
  {
    std::cout << "Initial seed candidate size is " << seeds.size() << std::endl;
  }
  std::sort(seeds.begin(), seeds.end(),
            [](const seed& a, const seed& b)
            { return a.ckeys.size() > b.ckeys.size(); });

  auto prunedSeeds = combineSeeds(seeds, clusterPositions);
  if (Verbosity() > 1)
  {
    std::cout << "Pruned seed size is " << prunedSeeds.size() << std::endl;
  }
  std::sort(prunedSeeds.begin(), prunedSeeds.end(),
            [](const seed& a, const seed& b)
            { return a.ckeys.size() > b.ckeys.size(); });

  auto finalSeeds = findIntersections(prunedSeeds);
  if (Verbosity() > 1)
  {
    std::cout << "final seeds are " << finalSeeds.size() << std::endl;
  }
  for (auto& seed1 : finalSeeds)
  {
    recalculateSeedLineParameters(seed1, clusterPositions, true);

    recalculateSeedLineParameters(seed1, clusterPositions, false);
  }
  auto chainedSeeds = chainSeeds(finalSeeds, clusterPositions);
  // auto chainedSeeds = finalSeeds;
  if (Verbosity() > 1)
  {
    std::cout << "Total seeds found is " << chainedSeeds.size() << std::endl;
  }
  int iseed = 0;

  auto longestseedit = chainedSeeds.begin();
  seed longestseed;
  if (longestseedit != chainedSeeds.end())
  {
    longestseed = chainedSeeds.front();
  }
  for (auto& seed_A : chainedSeeds)
  {
    // if mvtx only, we are interested only in seeds with > 3 hits
    if (m_trackerId == TrkrDefs::TrkrId::mvtxId && seed_A.ckeys.size() < 4)
    {
      continue;
    }
    if (m_analysis)
    {
      float seed_data[] = {
          (float) m_event,
          (float) iseed,
          (float) seed_A.ckeys.size(),
          seed_A.xyintercept,
          seed_A.xyslope,
          seed_A.xzintercept,
          seed_A.xzslope,
          longestseed.xyintercept,
          longestseed.xyslope,
          longestseed.xzintercept,
          longestseed.xzslope};
      m_tup->Fill(seed_data);
    }
    auto svtxseed = std::make_unique<TrackSeed_v2>();
    for (auto& key : seed_A.ckeys)
    {
      svtxseed->insert_cluster_key(key);
    }
    // if (m_trackerId == TrkrDefs::TrkrId::mvtxId){
    //   svtxseed->circleFitByTaubin(clusterPositions, 0, 3);
    //   svtxseed->lineFit(clusterPositions, 0, 3);
    // }
    m_seedContainer->insert(svtxseed.get());
    ++iseed;
  }
  if (Verbosity() > 1)
  {
    std::cout << "Final n seeds: " << m_seedContainer->size() << std::endl;
  }
  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}
PHCosmicSeeder::SeedVector PHCosmicSeeder::findIntersections(PHCosmicSeeder::SeedVector& initialSeeds)
{
  std::set<int> seedsToDelete;
  //! combine seeds with common ckeys
  for (unsigned int i = 0; i < initialSeeds.size(); ++i)
  {
    auto& seed1 = initialSeeds[i];
    for (unsigned int j = i; j < initialSeeds.size(); ++j)
    {
      if (i == j)
      {
        continue;
      }
      auto& seed2 = initialSeeds[j];
      std::vector<TrkrDefs::cluskey> intersection;
      unsigned int intersection_size_limit = 3;
      if (m_trackerId == TrkrDefs::TrkrId::mvtxId)
      {
        intersection_size_limit = 2;
      }
      std::set_intersection(seed1.ckeys.begin(), seed1.ckeys.end(),
                            seed2.ckeys.begin(), seed2.ckeys.end(), std::back_inserter(intersection));
      if (intersection.size() > intersection_size_limit)
      {
        //! If they share at least 4 (3 for MVTX only) clusters they are likely the same track,
        //! so merge and delete
        for (auto key : seed2.ckeys)
        {
          seed1.ckeys.insert(key);
        }
        seedsToDelete.insert(j);
      }
    }
  }

  SeedVector finalSeeds;
  unsigned int clus_per_seed_size = 4;
  if (m_trackerId == TrkrDefs::TrkrId::mvtxId)
  {
    clus_per_seed_size = 3;
  }
  for (unsigned int i = 0; i < initialSeeds.size(); ++i)
  {
    if (seedsToDelete.find(i) != seedsToDelete.end())
    {
      continue;
    }
    if (initialSeeds[i].ckeys.size() < clus_per_seed_size)
    {
      continue;
    }
    finalSeeds.push_back(initialSeeds[i]);
  }
  return finalSeeds;
}
PHCosmicSeeder::SeedVector PHCosmicSeeder::chainSeeds(PHCosmicSeeder::SeedVector& initialSeeds,
                                                      PositionMap& clusterPositions)
{
  PHCosmicSeeder::SeedVector returnseeds;
  std::set<int> seedsToDelete;
  if (Verbosity() > 1)
  {
    std::cout << "Chaining seeds, n seeds =  " << initialSeeds.size() << std::endl;
  }
  for (unsigned int i = 0; i < initialSeeds.size(); ++i)
  {
    auto& seed1 = initialSeeds[i];
    for (unsigned int j = i; j < initialSeeds.size(); ++j)
    {
      if (i == j)
      {
        continue;
      }
      auto& seed2 = initialSeeds[j];
      recalculateSeedLineParameters(seed1, clusterPositions, true);
      recalculateSeedLineParameters(seed2, clusterPositions, true);
      recalculateSeedLineParameters(seed1, clusterPositions, false);
      recalculateSeedLineParameters(seed2, clusterPositions, false);

      float longestxyslope = seed1.xyslope;
      float longestxyint = seed1.xyintercept;
      float longestxzslope = seed1.xzslope;
      float longestxzint = seed1.xzintercept;
      if (seed1.ckeys.size() < seed2.ckeys.size())
      {
        longestxyint = seed2.xyintercept;
        longestxyslope = seed2.xyslope;
        longestxzint = seed2.xzintercept;
        longestxzslope = seed2.xzslope;
      }

      float pdiff_tol = 1.0;
      if (m_trackerId == TrkrDefs::TrkrId::mvtxId)
      {
        pdiff_tol = 0.25;
      }
      float const pdiff = std::abs((seed1.xyslope - seed2.xyslope) / longestxyslope);
      float const pdiff2 = std::abs((seed1.xyintercept - seed2.xyintercept) / longestxyint);
      float const pdiff3 = std::abs((seed1.xzintercept - seed2.xzintercept) / longestxzint);
      float const pdiff4 = std::abs((seed1.xzslope - seed2.xzslope) / longestxzslope);
      if (Verbosity() > 1)
      {
        std::cout << "pdiff1,2,3,4 " << pdiff << ", " << pdiff2 << ", " << pdiff3 << ", " << pdiff4 << std::endl;
      }
      if (pdiff < pdiff_tol && pdiff2 < pdiff_tol && pdiff3 < pdiff_tol && pdiff4 < pdiff_tol)
      {
        seedsToDelete.insert(j);
        for (auto& key : seed2.ckeys)
        {
          seed1.ckeys.insert(key);
        }
      }
    }
  }
  for (unsigned int i = 0; i < initialSeeds.size(); ++i)
  {
    if (seedsToDelete.find(i) != seedsToDelete.end())
    {
      continue;
    }

    // if mvtx only, require each seed have at least one hit in each layer
    if (m_trackerId == TrkrDefs::TrkrId::mvtxId)
    {
      std::vector<bool> contains_layer = {false, false, false};
      for (auto& key : initialSeeds[i].ckeys)
      {
        contains_layer[TrkrDefs::getLayer(key)] = true;
      }
      if (std::find(contains_layer.begin(), contains_layer.end(), false) != contains_layer.end())
      {
        continue;
      }
    }

    returnseeds.push_back(initialSeeds[i]);
  }
  if (returnseeds.size() > 1 && m_trackerId == TrkrDefs::TrkrId::mvtxId)
  {
    PHCosmicSeeder::SeedVector emptyseeds;
    return emptyseeds;
  }

  return returnseeds;
}
PHCosmicSeeder::SeedVector PHCosmicSeeder::combineSeeds(PHCosmicSeeder::SeedVector& initialSeeds,
                                                        PHCosmicSeeder::PositionMap& clusterPositions)
{
  SeedVector prunedSeeds;
  std::set<int> seedsToDelete;

  for (unsigned int i = 0; i < initialSeeds.size(); ++i)
  {
    for (unsigned int j = i; j < initialSeeds.size(); ++j)
    {
      if (i == j)
      {
        continue;
      }
      auto& seed1 = initialSeeds[i];
      auto& seed2 = initialSeeds[j];

      // recalculate seed parameters
      recalculateSeedLineParameters(seed1, clusterPositions, true);
      recalculateSeedLineParameters(seed2, clusterPositions, true);
      recalculateSeedLineParameters(seed1, clusterPositions, false);
      recalculateSeedLineParameters(seed2, clusterPositions, false);

      float slope_tol = 0.1;
      float incept_tol = 3.0;
      if (m_trackerId == TrkrDefs::TrkrId::mvtxId)
      {
        slope_tol = 0.1;
        incept_tol = 0.5;
      }
      if (Verbosity() > 4)
      {
        std::cout << "xy slope diff: " << seed1.xyslope - seed2.xyslope << ", xy int diff   " << seed1.xyintercept - seed2.xyintercept << ", xz slope diff " << seed1.xzslope - seed2.xzslope << ", xz int diff   " << seed1.xzintercept - seed2.xzintercept << std::endl;
      }
      //! These values are tuned on the cosmic data
      if (std::abs(seed1.xyslope - seed2.xyslope) < slope_tol &&
          std::abs(seed1.xyintercept - seed2.xyintercept) < incept_tol &&
          std::abs(seed1.xzslope - seed2.xzslope) < slope_tol &&
          std::abs(seed1.xzintercept - seed2.xzintercept) < incept_tol)
      {
        for (auto& key : seed2.ckeys)
        {
          seed1.ckeys.insert(key);
        }
        seedsToDelete.insert(j);
      }
    }
  }
  if (Verbosity() > 4)
  {
    std::cout << "seeds to delete size is " << seedsToDelete.size() << std::endl;
  }
  for (unsigned int i = 0; i < initialSeeds.size(); ++i)
  {
    if (seedsToDelete.find(i) != seedsToDelete.end())
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
  //  int seednum = 0;
  double dist_check = 2.0;
  if (m_trackerId == TrkrDefs::TrkrId::mvtxId)
  {
    dist_check = 2.5;
  }
  for (auto& [key1, pos1] : clusterPositions)
  {
    for (auto& [key2, pos2] : clusterPositions)
    {
      if (key1 == key2)
      {
        continue;
      }
      // make a cut on clusters to at least be close to each other within a few cm
      float const dist = (pos2 - pos1).norm();
      if (m_trackerId == TrkrDefs::TrkrId::mvtxId && (TrkrDefs::getLayer(key1) == TrkrDefs::getLayer(key2)))
      {
        continue;
      }
      if (dist > dist_check)
      {
        continue;
      }
      if (Verbosity() > 5)
      {
        std::cout << "checking keys " << key1 << ", " << key2 << " with dist apart " << dist << std::endl;
        std::cout << "positions are " << pos1.transpose() << "    ,     "
                  << pos2.transpose() << std::endl;
      }
      PHCosmicSeeder::seed doub;

      doub.xyslope = (pos2.y() - pos1.y()) / (pos2.x() - pos1.x());
      doub.xyintercept = pos1.y() - doub.xyslope * pos1.x();
      doub.xzslope = (pos2.z() - pos1.z()) / (pos2.x() - pos1.x());
      doub.xzintercept = pos1.z() - pos1.x() * doub.xzslope;
      doub.yzslope = (pos2.z() - pos1.z()) / (pos2.y() - pos1.y());
      doub.yzintercept = pos1.z() - pos1.y() * doub.yzslope;

      keys.insert(key1);
      keys.insert(key2);
      doub.ckeys = keys;
      keys.clear();
      seeds.push_back(doub);
      //      seednum++;
    }
  }
  if (Verbosity() > 2)
  {
    std::cout << "doublet sizes " << seeds.size() << std::endl;
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
      float const dist1 = (pos1 - pos).norm();
      float const dist2 = (pos2 - pos).norm();
      float dist12_check = 2.;
      if (m_trackerId == TrkrDefs::TrkrId::mvtxId)
      {
        dist12_check = 3.5;
      }
      if (dist1 < dist2)
      {
        if (dist1 > dist12_check)
        {
          continue;
        }
      }
      else
      {
        if (dist2 > dist12_check)
        {
          continue;
        }
      }

      float const predy = dub.xyslope * pos.x() + dub.xyintercept;
      float const predz = dub.xzslope * pos.x() + dub.xzintercept;
      float const predz2 = dub.yzslope * pos.y() + dub.yzintercept;
      if (Verbosity() > 2)
      {
        std::cout << "testing ckey " << key << " with box dca "
                  << predy << ", " << pos.transpose() << " and " << predz << ", " << pos.z()
                  << std::endl;
      }
      if (fabs(predy - pos.y()) < m_xyTolerance)
      {
        if (m_trackerId == TrkrDefs::TrkrId::mvtxId && (fabs(predz - pos.z()) > 0.3 || fabs(predz2 - pos.z()) > 0.3))
        {
          continue;
        }
        if (Verbosity() > 2)
        {
          std::cout << "   adding ckey " << key << " with box dca "
                    << predy << ", " << pos.y() << " and " << predz << ", " << pos.z()
                    << std::endl;
          std::cout << "       with pos " << pos.x() << ", " << pos.y() << ", " << pos.z()
                    << std::endl;
        }
        dub.ckeys.insert(key);
      }
    }
  }

  /// erase doublets
  std::set<int> seedsToDelete;
  for (unsigned int i = 0; i < seeds.size(); ++i)
  {
    auto seed_A = seeds[i];
    if (Verbosity() > 2 && seed_A.ckeys.size() > 2)
    {
      std::cout << "keys in seed " << std::endl;
      for (auto& key : seed_A.ckeys)
      {
        std::cout << key << ", ";
      }
      std::cout << "seed xy slope: " << seed_A.xyslope << std::endl;
      std::cout << std::endl
                << "seed pos " << std::endl;
      for (auto& key : seed_A.ckeys)
      {
        std::cout << clusterPositions.find(key)->second.transpose() << std::endl;
      }
      std::cout << "Done printing seed info" << std::endl;
    }
    if (seed_A.ckeys.size() < 3)
    {
      seedsToDelete.insert(i);
    }
  }

  PHCosmicSeeder::SeedVector returnSeeds;
  for (unsigned int i = 0; i < seeds.size(); i++)
  {
    if (seedsToDelete.find(i) != seedsToDelete.end())
    {
      continue;
    }
    returnSeeds.push_back(seeds[i]);
  }
  return returnSeeds;
}
void PHCosmicSeeder::recalculateSeedLineParameters(seed& seed_A,
                                                   PHCosmicSeeder::PositionMap& clusters, bool isXY)
{
  float avgx = 0;
  float avgy = 0;
  PHCosmicSeeder::PositionMap seedClusters;
  for (auto& key : seed_A.ckeys)
  {
    auto glob = clusters.find(key)->second;
    if (isXY)
    {
      avgx += glob.x();
      avgy += glob.y();
    }
    else
    {
      avgx += glob.x();
      avgy += glob.z();
    }
    seedClusters.insert(std::make_pair(key, glob));
  }

  avgx /= seed_A.ckeys.size();
  avgy /= seed_A.ckeys.size();
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
      num += (glob.x() - avgx) * (glob.z() - avgy);
      denom += square(glob.x() - avgx);
    }
  }
  if (isXY)
  {
    seed_A.xyslope = num / denom;
    seed_A.xyintercept = avgy - seed_A.xyslope * avgx;
  }
  else
  {
    seed_A.xzslope = num / denom;
    seed_A.xzintercept = avgy - seed_A.xzslope * avgx;
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
int PHCosmicSeeder::End(PHCompositeNode* /*unused*/)
{
  if (m_outfile)
  {
    m_outfile->cd();
    m_tup->Write();
    m_outfile->Close();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}