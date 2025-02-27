#include "SingleInttPoolInput.h"
#include "intt_pool.h"

#include "Fun4AllStreamingInputManager.h"
#include "InputManagerType.h"

#include <ffarawobjects/InttRawHitContainerv2.h>
#include <ffarawobjects/InttRawHitv2.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>

#include <algorithm>  // for max
#include <cstdint>    // for uint64_t
#include <cstdlib>    // for exit
#include <iostream>   // for operator<<, basic_o...
#include <set>
#include <memory>
#include <utility>  // for pair

SingleInttPoolInput::SingleInttPoolInput(const std::string &name)
  : SingleStreamingInput(name)
{
  SubsystemEnum(InputManagerType::INTT);
  plist = new Packet *[1];
  m_rawHitContainerName = "INTTRAWHIT";
}

SingleInttPoolInput::~SingleInttPoolInput()
{
  delete[] plist;
  for (auto iter : poolmap)
  {
    if (Verbosity() > 2)
    {
      std::cout << "deleting intt pool for id  " << (iter.second)->getIdentifier() << std::endl;
    }
    delete (iter.second);
  }
}

void SingleInttPoolInput::FillPool(const uint64_t minBCO)
{
  if (AllDone())  // no more files and all events read
  {
    return;
  }
  while (GetEventiterator() == nullptr)  // at startup this is a null pointer
  {
    if (!OpenNextFile())
    {
      AllDone(1);
      return;
    }
  }

  //  std::set<uint64_t> saved_beamclocks;
  while (GetSomeMoreEvents(0))
  {
    Event *evt = GetEventiterator()->getNextEvent();
    while (!evt)
    {
      fileclose();
      if (!OpenNextFile())
      {
        AllDone(1);
        return;
      }
      evt = GetEventiterator()->getNextEvent();
    }
    if (Verbosity() > 2)
    {
      std::cout << "Fetching next Event" << evt->getEvtSequence() << std::endl;
    }
    RunNumber(evt->getRunNumber());
    if (GetVerbosity() > 1)
    {
      evt->identify();
    }
    // not interested in special events, really
    if (evt->getEvtType() != DATAEVENT)
    {
      m_NumSpecialEvents++;
      if (evt->getEvtType() == ENDRUNEVENT)
      {
        std::cout << "End run flag for INTT found, remaining INTT data is corrupted" << std::endl;
        delete evt;
        AllDone(1);
        return;
      }
      delete evt;
      continue;
    }

    int EventSequence = evt->getEvtSequence();
    int npackets = evt->getPacketList(plist, 1);

    if (npackets > 1)
    {
      exit(1);
    }
    if (m_SkipEarlyEvents)
    {
      for (int i = 0; i < npackets; i++)
      {
        int numBCOs = plist[i]->iValue(0, "NR_BCOS");
        for (int j = 0; j < numBCOs; j++)
        {
          uint64_t bco = plist[i]->lValue(j, "BCOLIST");
          if (bco < minBCO)
          {
            continue;
          }
          m_SkipEarlyEvents = false;
        }
      }
    }
    if (m_SkipEarlyEvents)
    {
      for (int i = 0; i < npackets; i++)
      {
        delete plist[i];
      }
      delete evt;
      continue;
    }
    for (int i = 0; i < npackets; i++)
    {
      if (Verbosity() > 2)
      {
        plist[i]->identify();
      }

      if (poolmap.find(plist[i]->getIdentifier()) == poolmap.end())  // we haven't seen this one yet
      {
        if (Verbosity() > 1)
        {
          std::cout << "starting new intt pool for packet " << plist[i]->getIdentifier() << std::endl;
        }
        poolmap[plist[i]->getIdentifier()] = new intt_pool(1000, 100);
        poolmap[plist[i]->getIdentifier()]->Verbosity(Verbosity());
        poolmap[plist[i]->getIdentifier()]->Name(std::to_string(plist[i]->getIdentifier()));
      }
      poolmap[plist[i]->getIdentifier()]->addPacket(plist[i]);

      delete plist[i];
    }

    delete evt;

    for (auto iter : poolmap)
    {
      intt_pool *pool = iter.second;  // less typing
      auto packet_id = pool->getIdentifier();
      if (pool->depth_ok())
      {
        int num_hits = pool->iValue(0, "NR_HITS");
        if (Verbosity() > 1)
        {
          std::cout << "Number of Hits: " << num_hits << " for packet "
                    << pool->getIdentifier() << std::endl;
        }

        int numBCOs = pool->iValue(0, "NR_BCOS");
        uint64_t largest_bco = 0;
        bool skipthis{true};

        for (int j = 0; j < numBCOs; j++)
        {
          uint64_t bco = pool->lValue(j, "BCOLIST");
          if (largest_bco < bco)
          {
            largest_bco = bco;
          }
          if (bco < minBCO)
          {
            continue;
          }
          skipthis = false;
          m_BclkStack.insert(bco);
          m_BclkStackPacketMap[packet_id].insert(bco);
        }

        int nFEEs = pool->iValue(0, "UNIQUE_FEES");
        for (int j = 0; j < nFEEs; j++)
        {
          int fee = pool->iValue(j, "FEE_ID");
          int nbcos = pool->iValue(fee, "FEE_BCOS");
          for (int k = 0; k < nbcos; k++)
          {
            uint64_t bco = pool->lValue(fee, k, "BCOVAL");
            if(bco < minBCO)
            {
              continue;
            }
            m_FeeGTML1BCOMap[fee].insert(bco);
          }
        }
        if (skipthis)
        {
          if (Verbosity() > 1)
          {
            std::cout << "largest bco: 0x" << std::hex << largest_bco << ", minbco 0x" << minBCO
                      << std::dec << ", evtno: " << EventSequence << std::endl;
          }
        }
        else
        {
          for (int j = 0; j < num_hits; j++)
          {
            uint64_t gtm_bco = pool->lValue(j, "BCO");
            if (gtm_bco < minBCO)
            {
              // std::cout << "dropping hit with bco 0x" << std::hex
              // 	      << gtm_bco << ", min bco: 0x" << minBCO
              // 	      << std::endl;
              continue;
            }
            auto newhit = std::make_unique<InttRawHitv2>();
            int FEE = pool->iValue(j, "FEE");
            newhit->set_packetid(pool->getIdentifier());
            newhit->set_fee(FEE);
            newhit->set_bco(gtm_bco);
            newhit->set_adc(pool->iValue(j, "ADC"));
            newhit->set_amplitude(pool->iValue(j, "AMPLITUDE"));
            newhit->set_chip_id(pool->iValue(j, "CHIP_ID"));
            newhit->set_channel_id(pool->iValue(j, "CHANNEL_ID"));
            newhit->set_word(pool->iValue(j, "DATAWORD"));
            newhit->set_FPHX_BCO(pool->iValue(j, "FPHX_BCO"));
            newhit->set_full_FPHX(pool->iValue(j, "FULL_FPHX"));
            newhit->set_full_ROC(pool->iValue(j, "FULL_ROC"));
            newhit->set_event_counter(pool->iValue(j, "EVENT_COUNTER"));
            gtm_bco += m_Rollover[FEE];

            if (gtm_bco < m_PreviousClock[FEE])
            {
              m_Rollover[FEE] += 0x10000000000;
              gtm_bco += 0x10000000000;  // rollover makes sure our bclks are ascending even if we roll over the 40 bit counter
            }
            m_PreviousClock[FEE] = gtm_bco;
            m_BeamClockFEE[gtm_bco].insert(FEE);
            m_FEEBclkMap[FEE] = gtm_bco;
            if (Verbosity() > 2)
            {
              std::cout << "evtno: " << EventSequence
                        << ", hits: " << j
                        << ", nr_hits: " << num_hits
                        << ", FEE: " << FEE
                        << ", bco: 0x" << std::hex << gtm_bco << std::dec
                        << ", min bco: 0x" << std::hex << minBCO << std::dec
                        << ", channel: " << newhit->get_channel_id()
                        << ", evt_counter: " << newhit->get_event_counter() << std::endl;
            }
            if (StreamingInputManager())
            {
              StreamingInputManager()->AddInttRawHit(gtm_bco, newhit.get());
            }
            m_InttRawHitMap[gtm_bco].push_back(newhit.release());
          }
        }
        //	    Print("FEEBCLK");
      }
      pool->next();
    }
  }
}

void SingleInttPoolInput::Print(const std::string &what) const
{
  if (what == "ALL" || what == "FEE")
  {
    for (const auto &bcliter : m_BeamClockFEE)
    {
      std::cout << "Beam clock 0x" << std::hex << bcliter.first << std::dec << std::endl;
      for (auto feeiter : bcliter.second)
      {
        std::cout << "FEM: " << feeiter << std::endl;
      }
    }
  }
  if (what == "ALL" || what == "FEEBCLK")
  {
    std::cout << "Printing last beamclock for every FEE" << std::endl;
    for (auto bcliter : m_FEEBclkMap)
    {
      std::cout << "FEE" << bcliter.first << " bclk: 0x"
                << std::hex << bcliter.second << std::dec << std::endl;
    }
  }
  if (what == "ALL" || what == "STORAGE")
  {
    for (const auto &bcliter : m_InttRawHitMap)
    {
      std::cout << "Beam clock 0x" << std::hex << bcliter.first << std::dec << std::endl;
      for (auto feeiter : bcliter.second)
      {
        std::cout << "fee: " << feeiter->get_fee()
                  << " at " << std::hex << feeiter << std::dec << std::endl;
      }
    }
  }
  if (what == "ALL" || what == "STACK")
  {
    for (auto &[packetid, bclkstack] : m_BclkStackPacketMap)
    {
      for (auto &bclk : bclkstack)
      {
        std::cout << "stacked bclk: 0x" << std::hex << bclk << std::dec << std::endl;
      }
    }
    for (auto iter : m_BclkStack)
    {
      std::cout << "stacked bclk: 0x" << std::hex << iter << std::dec << std::endl;
    }
  }
}

void SingleInttPoolInput::CleanupUsedPackets(const uint64_t bclk)
{
  m_BclkStack.erase(m_BclkStack.begin(), m_BclkStack.upper_bound(bclk));
  m_BeamClockFEE.erase(m_BeamClockFEE.begin(), m_BeamClockFEE.upper_bound(bclk));
  for(auto it = m_InttRawHitMap.begin(); it != m_InttRawHitMap.end() && (it->first <= bclk); it = m_InttRawHitMap.erase(it))
  {
    for( const auto& rawhit : it->second)
    {
      delete rawhit;
    }
  }
  m_InttRawHitMap.erase(m_InttRawHitMap.begin(), m_InttRawHitMap.upper_bound(bclk));
}

bool SingleInttPoolInput::CheckPoolDepth(const uint64_t bclk)
{
  for (auto iter : m_FEEBclkMap)
  {
    if (Verbosity() > 2)
    {
      std::cout << "my bclk 0x" << std::hex << iter.second
                << " req: 0x" << bclk << std::dec << std::endl;
    }
    if (iter.second < bclk)
    {
      if (Verbosity() > 1)
      {
        std::cout << "FEE " << iter.first << " beamclock 0x" << std::hex << iter.second
                  << " smaller than req bclk: 0x" << bclk << std::dec << std::endl;
      }
      return false;
    }
  }
  return true;
}

void SingleInttPoolInput::ClearCurrentEvent()
{
  // called interactively, to get rid of the current event
  uint64_t currentbclk = *(m_BclkStackPacketMap.begin()->second).begin();
  //  std::cout << "clearing bclk 0x" << std::hex << currentbclk << std::dec << std::endl;
  CleanupUsedPackets(currentbclk);
  // m_BclkStack.erase(currentbclk);
  // m_BeamClockFEE.erase(currentbclk);
  return;
}

bool SingleInttPoolInput::GetSomeMoreEvents(const uint64_t ibclk)
{
  if (AllDone())
  {
    return false;
  }
  if (poolmap.empty())
  {
    //      std::cout << "GetSomeMoreEvents poolmap empty, ret true" << std::endl;
    return true;
  }
  // for (auto iter : poolmap)
  // {
  //   if (!iter.second->depth_ok())
  //   {
  //   std::cout << "GetSomeMoreEvents depth not ok, ret true" << std::endl;
  // 	return true;
  //   }
  // }
  uint64_t localbclk = ibclk;
  if (ibclk == 0)
  {
    if (m_InttRawHitMap.empty())
    {
      //      std::cout << "GetSomeMoreEvents hitmap empty, ret true" << std::endl;
      return true;
    }
    localbclk = m_InttRawHitMap.begin()->first;
  }

  std::set<int> toerase;
  for (auto bcliter : m_FEEBclkMap)
  {
    if (bcliter.second <= localbclk)
    {
      uint64_t highest_bclk = m_InttRawHitMap.rbegin()->first;
      if ((highest_bclk - m_InttRawHitMap.begin()->first) < MaxBclkDiff())
      {
        // std::cout << "FEE " << bcliter.first << " bclk: "
        // 		<< std::hex << bcliter.second << ", req: " << localbclk
        // 		<< std::dec << std::endl;
        return true;
      }
      else
      {
        std::cout << PHWHERE << Name() << ": erasing FEE " << bcliter.first
                  << " with stuck bclk: 0x" << std::hex << bcliter.second
                  << " current bco range: 0x" << m_InttRawHitMap.begin()->first
                  << ", to: 0x" << highest_bclk << ", delta: " << std::dec
                  << (highest_bclk - m_InttRawHitMap.begin()->first)
                  << std::dec << std::endl;
        toerase.insert(bcliter.first);
      }
    }
  }
  for (auto iter : toerase)
  {
    m_FEEBclkMap.erase(iter);
  }
  //  std::cout << "GetSomeMoreEvents ret false" << std::endl;
  return false;
}

void SingleInttPoolInput::CreateDSTNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  PHNodeIterator iterDst(dstNode);
  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", "INTT"));
  if (!detNode)
  {
    detNode = new PHCompositeNode("INTT");
    dstNode->addNode(detNode);
  }
  InttRawHitContainer *intthitcont = findNode::getClass<InttRawHitContainer>(detNode, m_rawHitContainerName);
  if (!intthitcont)
  {
    intthitcont = new InttRawHitContainerv2();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(intthitcont, m_rawHitContainerName, "PHObject");
    detNode->addNode(newNode);
  }
}
//_______________________________________________________

void SingleInttPoolInput::ConfigureStreamingInputManager()
{
  if (StreamingInputManager())
  {
    StreamingInputManager()->SetInttBcoRange(m_BcoRange);
    StreamingInputManager()->SetInttNegativeBco(m_NegativeBco);
  }
}
