
#include "AzimuthalSeeder.h"

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

#include <TFile.h>
#include <TH2.h>

#include <cmath>
//____________________________________________________________________________..
AzimuthalSeeder::AzimuthalSeeder(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int AzimuthalSeeder::Init(PHCompositeNode * /*unused*/)
{
  file = std::make_unique<TFile>("ASoutfile.root", "recreate");
  h_phi = new TH2F("h_phi", "h_phi", 1000, -3.2, 3.2, 1000, -3.2, 3.2);
  h_phi2 = new TH2F("h_phi2", "h_phi2", 1000, -3.2, 3.2, 1000, -3.2, 3.2);
  h_phi3 = new TH2F("h_phi3", "h_phi3", 1000, -3.2, 3.2, 1000, -3.2, 3.2);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int AzimuthalSeeder::InitRun(PHCompositeNode *topNode)
{
  int ret = createNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }
  getNodes(topNode);
  return ret;
}

//____________________________________________________________________________..
int AzimuthalSeeder::process_event(PHCompositeNode * /*unused*/)
{
  AzimuthalSeeder::PositionMap clusterPositions[4];
  std::set<TrkrDefs::TrkrId> detectors;
  detectors.insert(TrkrDefs::TrkrId::mvtxId);
  for (const auto &det : detectors)
  {
    for (const auto &layer : {0, 1, 2})
    {
      for (const auto &hitsetkey : m_clusterContainer->getHitSetKeys(det, layer))
      {
        auto range = m_clusterContainer->getClusters(hitsetkey);
        for (auto citer = range.first; citer != range.second; ++citer)
        {
          const auto ckey = citer->first;
          const auto cluster = citer->second;
          const auto global = m_tGeometry->getGlobalPosition(ckey, cluster);
          clusterPositions[layer].insert(std::make_pair(ckey, global));
        }
      }
    }
  }

  SeedVector seeds;
  for (auto iter = clusterPositions[0].begin(); iter != clusterPositions[0].end(); ++iter)
  {
    auto key0 = iter->first;
    auto pos0 = iter->second;
    for (auto &iter2 : clusterPositions[1])
    {
      auto key1 = iter2.first;
      auto pos1 = iter2.second;
      if (key0 == key1)
      {
        continue;
      }
      float phi0 = atan2(pos0.y(), pos0.x());
      float phi1 = atan2(pos1.y(), pos1.x());

      h_phi->Fill(phi0, phi0 - phi1);
      if (std::fabs(phi0 - phi1) < 0.1)
      {
        seed s;
        s.ckeys.push_back(key0);
        s.ckeys.push_back(key1);
        s.globpos.push_back(pos0);
        s.globpos.push_back(pos1);
        seeds.push_back(s);
      }
    }
  }
  for (auto &s : seeds)
  {
    for (const auto &[key2, pos2] : clusterPositions[2])
    {
      float phi2 = atan2(pos2.y(), pos2.x());
      float phi1 = atan2(s.globpos[1].y(), s.globpos[1].x());
      float phi0 = atan2(s.globpos[0].y(), s.globpos[0].x());
      h_phi2->Fill(phi0, phi0 - phi2);
      h_phi3->Fill(phi2, phi1 - phi2);
      if (std::fabs(phi0 - phi2) < 0.1 && std::fabs(phi1 - phi2) < 0.1)
      {
        s.ckeys.push_back(key2);
        s.globpos.push_back(pos2);
      }
    }
  }

  SeedVector finalseeds;
  for (auto &s : seeds)
  {
    TrackFitUtils::position_vector_t xypoints, rzpoints;
    for (auto &pos : s.globpos)
    {
      xypoints.push_back(std::make_pair(pos.x(), pos.y()));
      float r = std::sqrt(pos.x() * pos.x() + pos.y() * pos.y());
      if (pos.y() < 0)
      {
        r = -r;
      }
      rzpoints.push_back(std::make_pair(pos.z(), r));
    }
    auto xyLineParams = TrackFitUtils::line_fit(xypoints);
    auto rzLineParams = TrackFitUtils::line_fit(rzpoints);
    float xySlope = std::get<0>(xyLineParams);
    float xyIntercept = std::get<1>(xyLineParams);
    float rzSlope = std::get<0>(rzLineParams);
    float rzIntercept = std::get<1>(rzLineParams);
    bool badseed = false;

    for (auto &pos : s.globpos)
    {
      float r = std::sqrt(pos.x() * pos.x() + pos.y() * pos.y());
      if (pos.y() < 0)
      {
        r = -r;
      }
      float distancexy = std::abs(xySlope * pos.x() - pos.y() + xyIntercept) / std::sqrt(xySlope * xySlope + 1);
      float distancerz = std::abs(rzSlope * pos.z() - r + rzIntercept) / std::sqrt(rzSlope * rzSlope + 1);
      if (distancexy > m_outlierLimit || distancerz > m_outlierLimit)
      {
        badseed = true;
        break;
      }
    }
    if (!badseed)
    {
      finalseeds.push_back(std::move(s));
    }
  }

  for (auto &map : {clusterPositions[1], clusterPositions[2]})
  {
    for (auto &[key, pos] : map)
    {
      clusterPositions[0].insert(std::make_pair(key, pos));
    }
  }
  if(Verbosity() > 2)
  {
    std::cout << "finalseed size " << finalseeds.size() << std::endl;
  }
  for (auto &s : finalseeds)
  {
    if (s.ckeys.size() < 3)
    {
      continue;
    }
    std::unique_ptr<TrackSeed_v2> si_seed = std::make_unique<TrackSeed_v2>();
    for (auto &key : s.ckeys)
    {
      si_seed->insert_cluster_key(key);
    }
    float avgphi = 0;
    for (auto &pos : s.globpos)
    {
      avgphi += atan2(pos.y(), pos.x());
    }
    avgphi /= s.globpos.size();
    si_seed->set_phi(avgphi);
    si_seed->circleFitByTaubin(clusterPositions[0], 0, 3);
    si_seed->lineFit(clusterPositions[0], 0, 3);
    if (Verbosity() > 3)
    {
      si_seed->identify();
    }
    m_seedContainer->insert(si_seed.get());
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int AzimuthalSeeder::End(PHCompositeNode * /*unused*/)
{
  if (m_outfile)
  {
    file->cd();
    h_phi->Write();
    h_phi2->Write();
    h_phi3->Write();
    file->Close();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int AzimuthalSeeder::createNodes(PHCompositeNode *topNode)
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

int AzimuthalSeeder::getNodes(PHCompositeNode *topNode)
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
