#include "SingleInttPoolInput.h"
#include "intt_pool.h"

#include "Fun4AllStreamingInputManager.h"

#include <ffarawobjects/InttRawHitContainerv1.h>
#include <ffarawobjects/InttRawHitv1.h>

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
#include <utility>  // for pair

SingleInttPoolInput::SingleInttPoolInput(const std::string &name)
  : SingleStreamingInput(name)
{
  SubsystemEnum(Fun4AllStreamingInputManager::INTT);
  plist = new Packet *[1];
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

void SingleInttPoolInput::FillPool(const unsigned int)
{
  if (AllDone())  // no more files and all events read
  {
    return;
  }
  while (GetEventiterator() == nullptr)  // at startup this is a null pointer
  {
    OpenNextFile();
  }

  std::set<uint64_t> saved_beamclocks;
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
      delete evt;
      continue;
    }

    int EventSequence = evt->getEvtSequence();
    int npackets = evt->getPacketList(plist, 1);

    if (npackets > 1)
    {
      exit(1);
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
      if (pool->depth_ok())
      {
        int num_hits = pool->iValue(0, "NR_HITS");
        if (Verbosity() > 1)
        {
          std::cout << "Number of Hits: " << num_hits << " for packet "
                    << pool->getIdentifier() << std::endl;
        }
        std::set<uint64_t> bclk_set;
        for (int j = 0; j < num_hits; j++)
        {
          InttRawHit *newhit = new InttRawHitv1();
          int FEE = pool->iValue(j, "FEE");
          uint64_t gtm_bco = pool->lValue(j, "BCO");
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

          gtm_bco += m_Rollover[FEE];
          bclk_set.insert(gtm_bco);
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
                      << ", channel: " << newhit->get_channel_id() << std::endl;
          }
          if (StreamingInputManager())
          {
            StreamingInputManager()->AddInttRawHit(gtm_bco, newhit);
          }
          m_InttRawHitMap[gtm_bco].push_back(newhit);
          m_BclkStack.insert(gtm_bco);
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
    for (auto iter : m_BclkStack)
    {
      std::cout << "stacked bclk: 0x" << std::hex << iter << std::dec << std::endl;
    }
  }
}

void SingleInttPoolInput::CleanupUsedPackets(const uint64_t bclk)
{
  std::vector<uint64_t> toclearbclk;
  for (const auto &iter : m_InttRawHitMap)
  {
    if (iter.first <= bclk)
    {
      for (auto pktiter : iter.second)
      {
        delete pktiter;
      }
      toclearbclk.push_back(iter.first);
    }
    else
    {
      break;
    }
  }
  for (auto iter : toclearbclk)
  {
    m_BclkStack.erase(iter);
    m_BeamClockFEE.erase(iter);
    m_InttRawHitMap.erase(iter);
  }
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
  uint64_t currentbclk = *m_BclkStack.begin();
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
    return true;
  }
  uint64_t localbclk = ibclk;
  if (ibclk == 0)
  {
    if (m_InttRawHitMap.empty())
    {
      return true;
    }
    localbclk = m_InttRawHitMap.begin()->first;
  }

  for (auto bcliter : m_FEEBclkMap)
  {
    if (bcliter.second <= localbclk)
    {
      // std::cout << "FEE " << bcliter.first << " bclk: "
      // 		<< std::hex << bcliter.second << ", req: " << localbclk
      // 		<< std::dec << std::endl;
      return true;
    }
  }
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
  InttRawHitContainer *intthitcont = findNode::getClass<InttRawHitContainer>(detNode, "INTTRAWHIT");
  if (!intthitcont)
  {
    intthitcont = new InttRawHitContainerv1();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(intthitcont, "INTTRAWHIT", "PHObject");
    detNode->addNode(newNode);
  }
}
