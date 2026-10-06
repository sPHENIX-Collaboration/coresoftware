#include "LaserAggregatedClusterizer.h"

#include "LaserEventInfo.h"

#include "LaserClusterHelper.h"

#include <trackbase/LaserAggregatedPad.h>
#include <trackbase/LaserAggregatedPadContainer.h>
#include <trackbase/LaserAggregatedPadContainerv1.h>
#include <trackbase/LaserAggregatedPadv1.h>
#include <trackbase/LaserCluster.h>
#include <trackbase/LaserClusterContainer.h>
#include <trackbase/LaserClusterContainerv1.h>
#include <trackbase/LaserClusterv5.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, getLayer

#include <ffaobjects/EventHeader.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco
#include <fun4all/Fun4AllUtils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>  // for PHIODataNode
#include <phool/PHNode.h>        // for PHNode
#include <phool/PHNodeIOManager.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TH2Poly.h>
#include <TMath.h>
#include <TFile.h>

#include <array>
#include <cmath>  // for sqrt, cos, sin
#include <format>
#include <iostream>
#include <fstream>
#include <limits>
#include <map>  // for _Rb_tree_cons...
#include <numeric>
#include <queue>
#include <set>
#include <string>
#include <tuple>
#include <utility>  // for pair
#include <vector>

namespace
{

  struct thread_data
  {
    uint8_t block = 0;
    std::map<padkey::key, double> adcMap;
    std::vector<LaserCluster *> cluster_vector;
    std::vector<TrkrDefs::cluskey> cluster_key_vector;
    int Verbosity = 0;
    pthread_mutex_t *mutex = nullptr;
  };

  void remove_hits(std::map<padkey::key, double> &adcMap, std::vector<padkey::key> &clusPads)
  {
    for (auto &pad : clusPads)
    {
      adcMap.erase(pad);
    }
  }

  void calc_cluster_parameter(std::vector<padkey::key>& clusPads, thread_data &my_data, uint8_t maxADCLayer)
  {
    
    unsigned int nHits = clusPads.size();
    if(nHits == 0)
    {
      return;
    }

    double layerSum = 0.0;
    double iphiSum = 0.0;

    double adcSum = 0.0;

    LaserCluster *clus = new LaserClusterv5;

    std::vector<uint8_t> usedLayer;
    std::vector<uint16_t> usedIPhi;

    float meanLayer = 0.0;
    float meanIPhi = 0.0;


    bool side = padkey::block_side(my_data.block);
    unsigned int sector = padkey::block_sector(my_data.block);

    for (auto &pad : clusPads)
    {
      uint8_t layer = padkey::get_layer(pad);
      uint16_t iphi = padkey::get_phibin(pad);
      double adc = my_data.adcMap.at(pad);

      bool foundLayer = false;
      for (int i : usedLayer)
      {
        if (layer == i)
        {
          foundLayer = true;
          break;
        }
      }

      if (!foundLayer)
      {
        usedLayer.push_back(layer);
      }

      bool foundIPhi = false;
      for (int i : usedIPhi)
      {
        if (iphi == i)
        {
          foundIPhi = true;
          break;
        }
      }

      if (!foundIPhi)
      {
        usedIPhi.push_back(iphi);
      }

      clus->addHitDouble(TpcDefs::genHitSetKey(layer, sector, side), TpcDefs::genHitKey(iphi, 30), adc);

      layerSum += 1.0 * layer * adc;
      iphiSum += 1.0 * iphi * adc;

      meanLayer += 1.0 * layer;
      meanIPhi += 1.0 * iphi;

      adcSum += adc;
    }

    if (nHits == 0 || clus->getNhitsDouble() == 0)
    {
      delete clus;
      return;
    }

    std::sort(usedLayer.begin(), usedLayer.end());
    std::sort(usedIPhi.begin(), usedIPhi.end());

    meanLayer = meanLayer / nHits;
    meanIPhi = meanIPhi / nHits;

    double sigmaLayer = 0.0;
    double sigmaIPhi = 0.0;

    double sigmaWeightedLayer = 0.0;
    double sigmaWeightedIPhi = 0.0;

    for (int i = 0; i < (int) clus->getNhitsDouble(); i++)
    {
      LaserClusterHitInfoDouble LCHI = clus->getHitDouble(i);
      uint8_t layer = TrkrDefs::getLayer(LCHI.hitsetkey);
      uint16_t iphi = TpcDefs::getPad(LCHI.hitkey);

      sigmaLayer += pow(layer - meanLayer, 2);
      sigmaIPhi += pow(iphi - meanIPhi, 2);

      sigmaWeightedLayer += LCHI.adc * pow(layer - (layerSum / adcSum), 2);
      sigmaWeightedIPhi += LCHI.adc * pow(iphi - (iphiSum / adcSum), 2);
    }

    clus->setNLayers(usedLayer.size());
    clus->setNIPhi(usedIPhi.size());
    clus->setNIT(0);
    clus->setSDLayer(sqrt(sigmaLayer / nHits));
    clus->setSDIPhi(sqrt(sigmaIPhi / nHits));
    clus->setSDIT(0.0);
    clus->setSDWeightedLayer(sqrt(sigmaWeightedLayer / adcSum));
    clus->setSDWeightedIPhi(sqrt(sigmaWeightedIPhi / adcSum));
    clus->setSDWeightedIT(0.0);

    const auto ckey = TrkrDefs::genClusKey(TpcDefs::genHitSetKey(maxADCLayer, sector, side), my_data.cluster_vector.size());    
    my_data.cluster_vector.push_back(clus);
    my_data.cluster_key_vector.push_back(ckey);

  }

  double clusterADCSum(std::map<padkey::key, double> &adcMap, const std::vector<padkey::key>& cluster)
  {
    double ADCSum = 0.0;
    for(const auto& pad : cluster)
    {
      ADCSum += adcMap.at(pad);
    }
    return ADCSum;
  }

  std::vector<padkey::key> trialFloodFill(std::map<padkey::key, double> &adcMap, padkey::key startingPad, int allowedAdjacentLayer)
  {
    std::vector<padkey::key> cluster;
    std::queue<padkey::key> queue;
    std::set<padkey::key> used;

    const uint8_t seedLayer = padkey::get_layer(startingPad);
    const uint16_t seedPhi = padkey::get_phibin(startingPad);

    queue.push(startingPad);
    used.insert(startingPad);

    while(!queue.empty())
    {
      auto current = queue.front();
      queue.pop();

      cluster.push_back(current);

      for(int dLayer = -1; dLayer <= 1; ++dLayer)
      {
        for(int dPhi = -1; dPhi <= 1; ++dPhi)
        {
          if(dLayer == 0 && dPhi == 0){ continue; }

          int nextLayer = padkey::get_layer(current) + dLayer;
          int nextPhi = padkey::get_phibin(current) + dPhi;

          if(nextLayer < padkey::kFirstTpcLayer || nextLayer >= padkey::kFirstTpcLayer + padkey::kNTpcLayers || nextPhi < 0){ continue; }

          padkey::key nextKey = padkey::genkey(nextLayer, padkey::get_side(current), padkey::get_sector(current), nextPhi);

          if(!adcMap.contains(nextKey)){ continue; }

          if(nextLayer != seedLayer && nextLayer != allowedAdjacentLayer){ continue; }

          if(std::abs(nextPhi - seedPhi) > 6){ continue; }

          if(used.contains(nextKey)){ continue; }

          queue.push(nextKey);
          used.insert(nextKey);

        }
      }
    }
    return cluster;
  }

  void ClusterModuleData(thread_data *my_data)
  {
    if (my_data->Verbosity > 2)
    {
      pthread_mutex_lock(my_data->mutex);
      std::cout << "clustering block: " << +my_data->block << "   side: " << +padkey::block_side(my_data->block) << "   sector: " << +padkey::block_sector(my_data->block) << "   module: " << +padkey::block_module(my_data->block) << std::endl;
      pthread_mutex_unlock(my_data->mutex);
    }

    while (!my_data->adcMap.empty())
    {
      auto iter = std::max_element(my_data->adcMap.begin(), my_data->adcMap.end(),
        [](const auto& pair1, const auto& pair2)
        {
          return pair1.second < pair2.second;
        }
      );

      padkey::key k = iter->first;
      uint8_t layer = padkey::get_layer(k);
      uint8_t module = padkey::get_module(k);

      if (my_data->Verbosity > 3)
      {
        pthread_mutex_lock(my_data->mutex);
        std::cout << "working on cluster " << my_data->cluster_vector.size() << "   side: " << +padkey::get_side(k) << "   sector: " << +padkey::get_sector(k) << "   module: " << +module << std::endl;
        pthread_mutex_unlock(my_data->mutex);
      }

      uint8_t minLayer;
      uint8_t maxLayer;

      if(module == 0)
      {
          minLayer = 7;
          maxLayer = 22;
      }
      else if(module == 1)
      {
          minLayer = 23;
          maxLayer = 38;
      }
      else
      {
          minLayer = 39;
          maxLayer = 54;
      }

      std::vector<padkey::key> lowerCluster;
      std::vector<padkey::key> upperCluster;
      
      if(layer > minLayer)
      {
        lowerCluster = trialFloodFill(my_data->adcMap, iter->first, layer-1);
      }
      if(layer < maxLayer)
      {
        upperCluster = trialFloodFill(my_data->adcMap, iter->first, layer+1);
      }
      
      std::vector<padkey::key> clusterPads;

      if(lowerCluster.empty() && !upperCluster.empty())
      {
        clusterPads = upperCluster;
      }
      else if(upperCluster.empty() && !lowerCluster.empty())
      {
        clusterPads = lowerCluster;
      }
      else
      {
      
        double lowerADCSum = clusterADCSum(my_data->adcMap, lowerCluster);
        double upperADCSum = clusterADCSum(my_data->adcMap, upperCluster);

        clusterPads = (lowerADCSum > upperADCSum) ? lowerCluster : upperCluster;
      }

      calc_cluster_parameter(clusterPads, *my_data, padkey::get_layer(k));

      remove_hits(my_data->adcMap, clusterPads);
    }
  }

  void *ClusterModule(void *threadarg)
  {
    auto *my_data = static_cast<thread_data *>(threadarg);
    ClusterModuleData(my_data);
    pthread_exit(nullptr);
  }
}  // namespace

LaserAggregatedClusterizer::LaserAggregatedClusterizer(const std::string &name)
  : SubsysReco(name)
{
}

int LaserAggregatedClusterizer::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Create the Cluster node if required
  auto *laserclusters = findNode::getClass<LaserClusterContainer>(dstNode, m_clusterNodeName);
  if (!laserclusters)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    laserclusters = new LaserClusterContainerv1;
    PHIODataNode<PHObject> *LaserClusterContainerNode =
        new PHIODataNode<PHObject>(laserclusters, m_clusterNodeName, "PHObject");
    DetNode->addNode(LaserClusterContainerNode);
  }

  m_clusterlist = findNode::getClass<LaserClusterContainer>(topNode, m_clusterNodeName);
  if (!m_clusterlist)
  {
    std::cout << PHWHERE << " ERROR: Can't find " << m_clusterNodeName << "." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  if (m_segfiles.empty())
  {
    std::cout << PHWHERE << " no segment files, call SetSegmentList or AddSegmentFile" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_aggregatedPads = new LaserAggregatedPadContainerv1;
  std::set<std::pair<int,int>> seen;
  int nMerged = 0;

  for (const auto &fname : m_segfiles)
  {
    PHNodeIOManager iman(fname, PHReadOnly, PHRunTree);
    if (!iman.isFunctional())
    {
      std::cout << PHWHERE << " no run tree in " << fname << std::endl;
      continue;
    }

    iman.selectObjectToRead("*",false);
    iman.selectObjectToRead("RUN#"+m_padContainerNodeName, true);

    PHCompositeNode *scratch = new PHCompositeNode("RUNSCRATCH");
    if (!iman.read(scratch))
    {
      delete scratch;
      continue;
    }

    if(Verbosity()>1)
    {
      std::cout << "working on segement " << nMerged << "/" << m_segfiles.size() << std::endl;
    }

    LaserAggregatedPadContainer *seg = findNode::getClass<LaserAggregatedPadContainer>(scratch, m_padContainerNodeName);
    if (seg)
    {
      std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(fname);
      if (seen.insert(runseg).second)
      {
        m_aggregatedPads->merge(seg);
        nMerged++;
      }
      else
      {
        std::cout << PHWHERE << " duplicate segment in list, skipped" << std::endl;
      }
    }
    else
    {
      std::cout << PHWHERE << m_padContainerNodeName << " missing in " << fname << std::endl;
    }
    delete scratch;

  }

  if (m_aggregatedPads->size() == 0)
  {
    delete m_aggregatedPads;
    m_aggregatedPads = nullptr;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  std::cout << "merged " << seen.size() << " segments, "
            << nMerged << " actually merged, "
            << m_aggregatedPads->size() << " pads, "
            << m_aggregatedPads->getNLaserEvents() << " laser events" << std::endl;

  auto *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << PHWHERE << "RUN Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  runNode->addNode(new PHIODataNode<PHObject>(m_aggregatedPads, std::format("{}Sum",m_padContainerNodeName), "PHObject"));

  return Fun4AllReturnCodes::EVENT_OK;
}

int LaserAggregatedClusterizer::process_event(PHCompositeNode* topNode)
{

  if(Verbosity()){ std::cout << "Already completed a round of clustering? " << m_done << std::endl; }

  if (m_done)
  {
    std::cout << PHWHERE << " clustering is single-shot, aborting run" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_done = true;

  if(!m_aggregatedPads || m_aggregatedPads->size() == 0)
  {
    std::cout << "LaserAggregatedClusterizer::process_event laser aggregated pads not found or were empty, aborting run" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  const double nEvents = m_aggregatedPads->getNLaserEvents();
  if(nEvents <= 0)
  {
    std::cout << "Zero laser events found in laser aggregated pads, aborting run" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  struct thread_pair_t
  {
    pthread_t thread{};
    thread_data data;
    bool started{false};
  };

  std::vector<thread_pair_t> threads;
  threads.reserve(padkey::kNBlocks);

  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  if (pthread_mutex_init(&m_threadlock, nullptr) != 0)
  {
    std::cout << std::endl
              << " mutex init failed" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  bool failed = false;
  for (uint8_t side = 0; side < padkey::kNSides; ++side)
  {
    for (uint8_t sec = 0; sec < padkey::kNSectors; ++sec)
    {
      for (uint8_t mod = 0; mod < padkey::kNModules; ++mod)
      {
        const uint8_t block = padkey::genblock(side, sec, mod);
        if (Verbosity() > 2)
        {
          std::cout << "making thread for block " << +block << std::endl;
          std::cout << "   side: " << +padkey::block_side(block) << "   sector: " << +padkey::block_sector(block) << "   module: " << +padkey::block_module(block) << std::endl;
        }

        thread_pair_t &thread_pair = threads.emplace_back();
        
        std::map<padkey::key, double> adcMap;
        auto range = m_aggregatedPads->getPadsInBlock(block);
        for(auto iter = range.first; iter != range.second; ++iter)
        {
          const double nHits = iter->second->getNHits();
          if(nHits / nEvents < m_nHitPerLaserEventMin){ continue; }
          const double ADC = iter->second->getAdc() / nEvents;
          const padkey::key k = iter->first;
          adcMap[k] = ADC;
        }

        std::vector<LaserCluster *> cluster_vector;
        std::vector<TrkrDefs::cluskey> cluster_key_vector;

        thread_pair.data.block = block;
        thread_pair.data.adcMap = adcMap;
        thread_pair.data.cluster_vector = cluster_vector;
        thread_pair.data.cluster_key_vector = cluster_key_vector;
        thread_pair.data.Verbosity = Verbosity();
        thread_pair.data.mutex = &m_threadlock;

        int rc;
        rc = pthread_create(&thread_pair.thread, &attr, ClusterModule, (void *) &thread_pair.data);

        if (rc != 0)
        {
          std::cerr << "Error:unable to create thread," << rc << std::endl;
          failed = true;
          continue;
        }
        thread_pair.started = true;
      }
    }
  }

 
  pthread_attr_destroy(&attr);

  LaserClusterHelper lch;
  lch.set_useDouble(true);
  lch.loadNodes(topNode); 
  lch.set_useZ(false);
  lch.set_useGlobal(true);

  PHG4TpcGeomContainer *geom_container = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  TFile *QAFile = nullptr;
  if(!m_QAName.empty())
  {
    QAFile = new TFile(m_QAName.c_str(),"RECREATE");
  }

  TH2Poly *h[2] = {nullptr, nullptr};
  if(QAFile)
  {
    for(int s=0; s<2; s++)
    {
      h[s] = new TH2Poly(std::format("hADC_perLaserFlash_{}",s ? "North" : "South").c_str(), std::format("ADC Per Laser Flash (nHits_perLaserFlash>=1) {}",s ? "North" : "South").c_str(),
                          100, 0.0, 2.0*TMath::Pi(), 100, 28, 78);
      h[s]->SetFloat(false);
      h[s]->SetDirectory(nullptr);
    }
  }


  for (const auto &thread_pair : threads)
  {
    if(!thread_pair.started){ continue; }
    int rc2 = pthread_join(thread_pair.thread, nullptr);
    if (rc2 != 0)
    {
      std::cout << "Error:unable to join," << rc2 << std::endl;
      failed = true;
      continue;
    }

    for (int index = 0; index < (int) thread_pair.data.cluster_vector.size(); ++index)
    {
      auto *cluster = thread_pair.data.cluster_vector[index];
      const auto ckey = thread_pair.data.cluster_key_vector[index];

      if(Verbosity() > 4){ cluster->identify(); }

      m_clusterlist->addClusterSpecifyKey(ckey, cluster);

      if(!m_QAName.empty() && geom_container)
      {
        if(Verbosity() > 3){ std::cout << "   working on cluster " << m_clusterlist->size() - 1 << std::endl; }
        int side = TpcDefs::getSide(ckey);
        for(int i=0; i<(int)cluster->getNhitsDouble(); i++)
        {
          if(Verbosity() > 4){ std::cout << "      working on hit " << i << std::endl; }
          LaserClusterHitInfoDouble LCHI = cluster->getHitDouble(i);
          int layer = TrkrDefs::getLayer(LCHI.hitsetkey);
          PHG4TpcGeom *layer_geom = geom_container->GetLayerCellGeom(layer);
          if(!layer_geom){ continue; }

          Acts::Vector3 global = lch.getHitPosition(LCHI.hitsetkey, LCHI.hitkey);

          double phiWidth = layer_geom->get_phistep() / 2.0;
          double layerWidth = layer_geom->get_thickness() / 2.0;

          double hitPhi = atan2(global[1], global[0]);
          while(hitPhi < 0){ hitPhi += 2*TMath::Pi(); }
          while(hitPhi > 2*TMath::Pi()){ hitPhi -= 2*TMath::Pi(); }
          double hitR = sqrt(global[0]*global[0] + global[1]*global[1]);      

          double phis[5] = {hitPhi - phiWidth,hitPhi + phiWidth,hitPhi + phiWidth,hitPhi - phiWidth,hitPhi - phiWidth};
          double rs[5] = {hitR - layerWidth, hitR - layerWidth, hitR + layerWidth, hitR + layerWidth, hitR - layerWidth};

          int bin = h[side]->AddBin(5, phis, rs);
          h[side]->SetBinContent(bin, LCHI.adc);
          if(Verbosity() > 3){ std::cout << "         added bin " << bin << " to TH2Poly " << (side ? "North" : "South") << " with content " << LCHI.adc << std::endl; }
        }
      }

    }
  }
  pthread_mutex_destroy(&m_threadlock);

  if(failed)
  {
    for(auto &hist : h)
    {
      delete hist;
    }
    if(QAFile){ QAFile->Close(); }
    delete QAFile;
    return Fun4AllReturnCodes::ABORTRUN;
  }


  if(!m_QAName.empty())
  {
    QAFile->cd();
    for(auto &hist : h)
    {
      hist->Write();
      delete hist;
    }
    QAFile->Close();
  }
  delete QAFile;

  
  return Fun4AllReturnCodes::EVENT_OK;
}

void LaserAggregatedClusterizer::SetSegmentList(const std::string &listfile)
{
  std::ifstream infile(listfile);
  if (!infile.is_open())
  {
    std::cout << PHWHERE << " could not open segment list " << listfile << std::endl;
    return;
  }

  std::string line;
  int nadded = 0;
  while (std::getline(infile, line))
  {
    // strip comments and surrounding whitespace
    const size_t hash = line.find('#');
    if (hash != std::string::npos){ line.erase(hash); }

    const size_t first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos){ continue; }          // blank or comment-only
    const size_t last = line.find_last_not_of(" \t\r\n");
    line = line.substr(first, last - first + 1);

    m_segfiles.push_back(line);
    ++nadded;
  }

  if (Verbosity() > 0)
  {
    std::cout << "LaserAggregatedClusterizer::SetSegmentList read " << nadded
              << " files from " << listfile << std::endl;
  }
}