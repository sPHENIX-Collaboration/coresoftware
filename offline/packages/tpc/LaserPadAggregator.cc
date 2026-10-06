#include "LaserPadAggregator.h"

#include "LaserEventInfo.h"

#include <trackbase/LaserAggregatedPad.h>
#include <trackbase/LaserAggregatedPadContainer.h>
#include <trackbase/LaserAggregatedPadContainerv1.h>
#include <trackbase/LaserAggregatedPadv1.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, getLayer
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <ffaobjects/EventHeader.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>  // for PHIODataNode
#include <phool/PHNode.h>        // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <array>
#include <cmath>  // for sqrt, cos, sin
#include <format>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>  // for _Rb_tree_cons...
#include <numeric>
#include <queue>
#include <set>
#include <string>
#include <tuple>
#include <utility>  // for pair
#include <vector>

#include <pthread.h>

namespace
{
  struct thread_data
  {
    std::vector<TrkrHitSet *> hitsets;
    std::vector<unsigned int> layers;
    bool side = false;
    unsigned int sector = 0;
    unsigned int module = 0;
    double adc_threshold = 74.4;
    int peakTimeBin = 325;
    int eventNum = 0;
    int Verbosity = 0;
    std::map<std::pair<uint8_t, uint16_t>, double> adcMap;
    std::map<std::pair<uint8_t, uint16_t>, double> hitMap;
  };

  pthread_mutex_t mythreadlock;

  void ProcessModuleData(thread_data *my_data)
  {
    if (my_data->Verbosity > 2)
    {
      pthread_mutex_lock(&mythreadlock);
      std::cout << "working on side: " << my_data->side << "   sector: " << my_data->sector << "   module: " << my_data->module << std::endl;
      pthread_mutex_unlock(&mythreadlock);
    }

    if (my_data->hitsets.empty())
    {
      return;
    }
    for (int i = 0; i < (int) my_data->hitsets.size(); i++)
    {
      auto *hitset = my_data->hitsets[i];
      unsigned int layer = my_data->layers[i];

      //TrkrDefs::hitsetkey hitsetkey = TpcDefs::genHitSetKey(layer, sector, (int) side);

      TrkrHitSet::ConstRange hitrangei = hitset->getHits();

      for (TrkrHitSet::ConstIterator hitr = hitrangei.first; hitr != hitrangei.second; ++hitr)
      {
        double_t fadc = hitr->second->getAdc();
        unsigned short adc = 0;
        if (fadc > my_data->adc_threshold)
        {
          adc = (unsigned short) fadc;
        }
        else
        {
          continue;
        }

        int iphi = TpcDefs::getPad(hitr->first);
        int it = TpcDefs::getTBin(hitr->first);

        if (std::abs(it - my_data->peakTimeBin) > 5)
        {
          continue;
        }

        my_data->adcMap[{layer, iphi}] += adc;
        my_data->hitMap[{layer, iphi}]++;
      }
    }
  } 

  void *ProcessModule(void *threadarg)
  {
    auto *my_data = static_cast<thread_data *>(threadarg);
    ProcessModuleData(my_data);
    pthread_exit(nullptr);
  }

} // namespace

LaserPadAggregator::LaserPadAggregator(const std::string &name)
  : SubsysReco(name)
{
}

int LaserPadAggregator::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << PHWHERE << "RUN Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Create the Cluster node if required
  auto *padlist = findNode::getClass<LaserAggregatedPadContainer>(runNode, m_padContainerNodeName);
  if (!padlist)
  {
    padlist = new LaserAggregatedPadContainerv1;
    PHIODataNode<PHObject> *LaserAggregatedPadContainerNode =
        new PHIODataNode<PHObject>(padlist, m_padContainerNodeName, "PHObject");
    runNode->addNode(LaserAggregatedPadContainerNode);
  }

  for(int s=0; s<2; s++)
  {
    for(int sec=0; sec<12; sec++)
    {
      for(int mod=0; mod<3; mod++)
      {
        m_adcMap[s][sec][mod].clear();
        m_hitMap[s][sec][mod].clear();
      }
    }
  }
  nLaserEvents = 0;
  

  return Fun4AllReturnCodes::EVENT_OK;
}

int LaserPadAggregator::process_event(PHCompositeNode *topNode)
{
  eventHeader = findNode::getClass<EventHeader>(topNode, "EventHeader");
  if (!eventHeader)
  {
    std::cout << PHWHERE << " EventHeader Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_event = eventHeader->get_EvtSequence();

  if (Verbosity() > 1)
  {
    std::cout << "LaserPadAggregator::process_event working on event " << m_event << std::endl;
  }

  m_laserEventInfo = findNode::getClass<LaserEventInfo>(topNode, "LaserEventInfo");
  if (!m_laserEventInfo)
  {
    std::cout << PHWHERE << "ERROR: Can't find node LaserEventInfo" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if ((eventHeader->get_RunNumber() > 66153 && !m_laserEventInfo->isGl1LaserEvent()) || (eventHeader->get_RunNumber() <= 66153 && !m_laserEventInfo->isLaserEvent()))
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if (Verbosity() > 1)
  {
    std::cout << "LaserPadAggregator::process_event laser event found" << std::endl;
  }

  // get node containing the digitized hits
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    std::cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  TrkrHitSetContainer::ConstRange hitsetrange = m_hits->getHitSets(TrkrDefs::TrkrId::tpcId);

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

  if (pthread_mutex_init(&mythreadlock, nullptr) != 0)
  {
    std::cout << std::endl
              << " mutex init failed" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  bool failed = false;
  for (unsigned int sec = 0; sec < 12; sec++)
  {
    for (int s = 0; s < 2; s++)
    {
      for (unsigned int mod = 0; mod < 3; mod++)
      {
        if (Verbosity() > 2)
        {
          std::cout << "making thread for side: " << s << "   sector: " << sec << "   module: " << mod << std::endl;
        }

        thread_pair_t &thread_pair = threads.emplace_back();

        std::vector<TrkrHitSet *> hitsets;
        std::vector<unsigned int> layers;

        for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
             hitsetitr != hitsetrange.second;
             ++hitsetitr)
        {
          unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
          int side = TpcDefs::getSide(hitsetitr->first);
          unsigned int sector = TpcDefs::getSectorId(hitsetitr->first);
          if (sector != sec || side != s)
          {
            continue;
          }
          if ((mod == 0 && (layer < 7 || layer > 22)) || (mod == 1 && (layer <= 22 || layer > 38)) || (mod == 2 && (layer <= 38 || layer > 54)))
          {
            continue;
          }

          TrkrHitSet *hitset = hitsetitr->second;

          hitsets.push_back(hitset);
          layers.push_back(layer);
        }

        thread_pair.data.hitsets = hitsets;
        thread_pair.data.layers = layers;
        thread_pair.data.side = (bool) s;
        thread_pair.data.sector = sec;
        thread_pair.data.module = mod;
        thread_pair.data.adc_threshold = m_adc_threshold;
        thread_pair.data.peakTimeBin = m_laserEventInfo->getPeakSample(s);
        thread_pair.data.eventNum = m_event;
        thread_pair.data.Verbosity = Verbosity();

        int rc;
        rc = pthread_create(&thread_pair.thread, &attr, ProcessModule, (void *) &thread_pair.data);

        if (rc != 0)
        {
          std::cout << "Error:unable to create thread," << rc << std::endl;
          failed = true;
          continue;
        }
        thread_pair.started = true;
      }
    }
  }

  pthread_attr_destroy(&attr);

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

    const auto& data = thread_pair.data;

    for(const auto& [pad, adc] : data.adcMap)
    {
      m_adcMap[data.side][data.sector][data.module][pad] += adc;
    }
    for(const auto& [pad, hits] : data.hitMap)
    {
      m_hitMap[data.side][data.sector][data.module][pad] += hits;
    }

  }


  threads.clear();
  pthread_mutex_destroy(&mythreadlock);

  if (failed) { return Fun4AllReturnCodes::ABORTRUN; }


  if (Verbosity() > 1)
  {
    std::cout << "LaserPadAggregator::process_event pads per sector: " << std::endl;
    std::cout << std::setw(5) << "side" << std::setw(7) << "sector" << std::setw(7) << "R1" << std::setw(7) << "R2" << std::setw(7) << "R3" << std::endl;
    for(int s=0; s<2; s++)
    {
      for(int sec=0; sec<12; sec++)
      {
        std::cout << std::setw(5) << (s==1 ? "North" : "South") << std::setw(7) << sec << std::setw(7) << m_adcMap[s][sec][0].size() << std::setw(7) << m_adcMap[s][sec][1].size() << std::setw(7) << m_adcMap[s][sec][2].size() << std::endl;
      }
    }
  }

  nLaserEvents++;

  return Fun4AllReturnCodes::EVENT_OK;
}

int LaserPadAggregator::End(PHCompositeNode *topNode)
{
  // get node for clusters
  m_padlist = findNode::getClass<LaserAggregatedPadContainer>(topNode, m_padContainerNodeName);
  if (!m_padlist)
  {
    std::cout << PHWHERE << " ERROR: Can't find " << m_padContainerNodeName << "." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_padlist->setNLaserEvents(nLaserEvents);

  for (unsigned int sec = 0; sec < 12; sec++)
  {
    for (int s = 0; s < 2; s++)
    {
      for (unsigned int mod = 0; mod < 3; mod++)
      {
        for(const auto& [pad, adc] : m_adcMap[s][sec][mod])
        {
          auto hitIt = m_hitMap[s][sec][mod].find(pad);
          if(hitIt != m_hitMap[s][sec][mod].end())
          {
            padkey::key new_padkey = padkey::genkey(pad.first, s, sec, pad.second);
            LaserAggregatedPad *newPad = new LaserAggregatedPadv1();
            newPad->addAdc(adc);
            newPad->addNHits(hitIt->second);
            m_padlist->addPad(new_padkey, newPad);
          }
        }
      }
    }
  }

  if (Verbosity() > 1)
  {
    std::cout << "LaserPadAggregator::End pads per sector: " << std::endl;
    std::cout << std::setw(5) << "side" << std::setw(7) << "sector" << std::setw(7) << "R1" << std::setw(7) << "R2" << std::setw(7) << "R3" << std::endl;
    for(int s=0; s<2; s++)
    {
      for(int sec=0; sec<12; sec++)
      {
        std::cout << std::setw(5) << (s==1 ? "North" : "South") << std::setw(7) << sec << std::setw(7) << m_adcMap[s][sec][0].size() << std::setw(7) << m_adcMap[s][sec][1].size() << std::setw(7) << m_adcMap[s][sec][2].size() << std::endl;
      }
    }
  }


  return Fun4AllReturnCodes::EVENT_OK;

}
