#include "SingleGl1PoolInputv2.h"

#include <fun4allraw/Fun4AllStreamingInputManager.h>
#include <fun4allraw/InputManagerType.h>
#include "Fun4AllStreamingLumiCountingInputManager.h"

#include <ffarawobjects/Gl1Packetv2.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/packet.h>  // for Packet

#include <cstdint>   // for uint64_t
#include <iostream>  // for operator<<, basic_ostream<...
#include <iterator>  // for reverse_iterator
#include <limits>    // for numeric_limits
#include <memory>
#include <set>
#include <utility>  // for pair

SingleGl1PoolInputv2::SingleGl1PoolInputv2(const std::string &name)
  : SingleStreamingInputv2(name)
{
  SubsystemEnum(InputManagerType::GL1);
}

SingleGl1PoolInputv2::~SingleGl1PoolInputv2()
{
  CleanupUsedPackets(std::numeric_limits<uint64_t>::max());
}

void SingleGl1PoolInputv2::FillPool(const unsigned int /*nbclks*/)
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
  while (GetSomeMoreEvents())
  {
    std::unique_ptr<Event> evt(GetEventiterator()->getNextEvent());
    while (!evt)
    {
      fileclose();
      if (!OpenNextFile())
      {
        AllDone(1);
        return;
      }
      evt.reset(GetEventiterator()->getNextEvent());
    }
    if (Verbosity() > 2)
    {
      std::cout << PHWHERE << "Fetching next Event" << evt->getEvtSequence() << std::endl;
    }
    if ((m_total_event == 0 && evt->getEvtType() == ENDRUNEVENT) ||
        (m_total_event != 0 && evt->getEvtSequence() - 2 == m_total_event))
    {
      m_alldone_flag = true;
      m_lastevent_flag = true;
    }
    if (evt->getEvtSequence() % 5000 == 0)
    {
      m_alldone_flag = true;
      m_lastevent_flag = true;
    }
    if (Verbosity() > 2)
    {
      if (m_alldone_flag)
      {
        std::cout << "gl1 all done is true" << std::endl;
      }
      // else{std::cout<<"gl1 all done is false"<<std::endl;}
    }
    RunNumber(evt->getRunNumber());
    if (GetVerbosity() > 1)
    {
      evt->identify();
    }
    if (evt->getEvtType() != DATAEVENT)
    {
      m_NumSpecialEvents++;
      if (evt->getEvtType() == ENDRUNEVENT)
      {
        AllDone(1);
        std::unique_ptr<Event> nextevt(GetEventiterator()->getNextEvent());
        if (nextevt)
        {
          std::cout << PHWHERE << " Found event after End Run Event " << std::endl;
          std::cout << "End Run Event identify: " << std::endl;
          evt->identify();
          std::cout << "Next event identify: " << std::endl;
          nextevt->identify();
        }
        return;
      }
      continue;
    }
    int EventSequence = evt->getEvtSequence();
    Packet *packet = evt->getPacket(14001);
    if (!packet)
    {
      std::cout << PHWHERE << "Packet 14001 is null ptr" << std::endl;
      evt->identify();
      m_alldone_flag = true;
      m_lastevent_flag = true;
      if (StreamingLumiInputManager())
      {
        StreamingLumiInputManager()->SetEndofEvent(m_alldone_flag, m_lastevent_flag);
        StreamingLumiInputManager()->SetEventNumber(EventSequence);
      }
      m_alldone_flag = false;
      m_lastevent_flag = false;
      continue;
    }
    if (Verbosity() > 1)
    {
      packet->identify();
    }

    Gl1Packet *newhit = new Gl1Packetv2();
    uint64_t gtm_bco = packet->lValue(0, "BCO");
    uint64_t bco_trim = gtm_bco & 0xFFFFFFFFFFU;
    //	std::cout<<bco_trim<<std::endl;
    m_BCOWindows[bco_trim] = std::make_pair(bco_trim - m_negative_bco_window, bco_trim + m_positive_bco_window);
    //    std::cout<<"BCO "<< m_BCOWindows.begin()->first<<" left "<<m_BCOWindows.begin()->second.first<<" right "<< m_BCOWindows.begin()->second.second<<std::endl;
    m_BCOBunchNumber[bco_trim] = packet->lValue(0, "BunchNumber");
    //  std::cout<<"BCO "<<bco_trim << " Bunch Number "<< m_BCOBunchNumber[bco_trim]<<std::endl;
    if (StreamingLumiInputManager())
    {
      StreamingLumiInputManager()->AddGl1Window(bco_trim, m_negative_bco_window, m_positive_bco_window);
      StreamingLumiInputManager()->AddGl1BunchNumber(bco_trim, m_BCOBunchNumber[bco_trim]);
      StreamingLumiInputManager()->SetEndofEvent(m_alldone_flag, m_lastevent_flag);
      StreamingLumiInputManager()->SetEventNumber(EventSequence);
      StreamingLumiInputManager()->SetNegativeWindow(m_negative_bco_window);
      StreamingLumiInputManager()->SetPositiveWindow(m_positive_bco_window);
    }
    if (evt->getEvtSequence() % 5000 == 0)
    {
      m_alldone_flag = false;
      m_lastevent_flag = false;
    }

    m_FEEBclkMap.insert(gtm_bco);
    newhit->setBCO(packet->lValue(0, "BCO"));
    newhit->setHitFormat(packet->getHitFormat());
    newhit->setIdentifier(packet->getIdentifier());
    newhit->setEvtSequence(EventSequence);
    newhit->setPacketNumber(packet->iValue(0));
    newhit->setBunchNumber(packet->lValue(0, "BunchNumber"));
    newhit->setTriggerInput(packet->lValue(0, "TriggerInput"));
    newhit->setLiveVector(packet->lValue(0, "LiveVector"));
    newhit->setScaledVector(packet->lValue(0, "ScaledVector"));
    newhit->setGTMBusyVector(packet->lValue(0, "GTMBusyVector"));
    for (int i = 0; i < 64; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        newhit->setScaler(i, j, packet->lValue(i, j));
      }
    }
    for (int i = 0; i < 12; i++)
    {
      newhit->setGl1pScaler(i, 0, packet->lValue(i, "GL1PRAW"));
      newhit->setGl1pScaler(i, 1, packet->lValue(i, "GL1PLIVE"));
      newhit->setGl1pScaler(i, 2, packet->lValue(i, "GL1PSCALED"));
    }
    if (Verbosity() > 2)
    {
      std::cout << PHWHERE << " Packet: " << packet->getIdentifier()
                << " evtno: " << EventSequence
                << ", bco: 0x" << std::hex << gtm_bco << std::dec
                << ", bunch no: " << packet->lValue(0, "BunchNumber")
                << std::endl;
      std::cout << PHWHERE << " RB Packet: " << newhit->getIdentifier()
                << " evtno: " << newhit->getEvtSequence()
                << ", bco: 0x" << std::hex << newhit->getBCO() << std::dec
                << ", bunch no: " << +newhit->getBunchNumber()
                << std::endl;
    }
    if (Verbosity() > 2)
    {
      std::cout << PHWHERE << "evtno: " << EventSequence
                << ", bco: 0x" << std::hex << gtm_bco << std::dec
                << std::endl;
    }
    if (StreamingLumiInputManager())
    {
      StreamingLumiInputManager()->AddGl1RawHit(gtm_bco, newhit);
    }

    m_Gl1RawHitMap[gtm_bco].push_back(newhit);
    m_BclkStack.insert(gtm_bco);

    delete packet;
  }
}

void SingleGl1PoolInputv2::Print(const std::string &what) const
{
  if (what == "ALL" || what == "FEEBCLK")
  {
    for (auto bcliter : m_FEEBclkMap)
    {
      std::cout << PHWHERE << " bclk: 0x"
                << std::hex << bcliter << std::dec << std::endl;
    }
  }
  if (what == "ALL" || what == "STORAGE")
  {
    for (const auto &bcliter : m_Gl1RawHitMap)
    {
      std::cout << PHWHERE << "Beam clock 0x" << std::hex << bcliter.first << std::dec << std::endl;
      for (auto feeiter : bcliter.second)
      {
        std::cout << PHWHERE << "fee: " << feeiter->getBCO()
                  << " at " << std::hex << feeiter << std::dec << std::endl;
      }
    }
  }
  if (what == "ALL" || what == "STACK")
  {
    for (auto iter : m_BclkStack)
    {
      std::cout << PHWHERE << "stacked bclk: 0x" << std::hex << iter << std::dec << std::endl;
    }
  }
}

void SingleGl1PoolInputv2::CleanupUsedPackets(const uint64_t bclk)
{
  std::vector<uint64_t> toclearbclk;
  for (const auto &iter : m_Gl1RawHitMap)
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
    m_FEEBclkMap.erase(iter);
    m_BclkStack.erase(iter);
    m_Gl1RawHitMap.erase(iter);
  }
}

bool SingleGl1PoolInputv2::CheckPoolDepth(const uint64_t bclk)
{
  // if (m_FEEBclkMap.size() < 10)
  // {
  //   std::cout << PHWHERE << "not all FEEs in map: " << m_FEEBclkMap.size() << std::endl;
  //   return true;
  // }
  for (auto iter : m_FEEBclkMap)
  {
    if (Verbosity() > 2)
    {
      std::cout << PHWHERE << "my bclk 0x" << std::hex << iter
                << " req: 0x" << bclk << std::dec << std::endl;
    }
    if (iter < bclk)
    {
      if (Verbosity() > 1)
      {
        std::cout << PHWHERE << "FEE " << iter << " beamclock 0x" << std::hex << iter
                  << " smaller than req bclk: 0x" << bclk << std::dec << std::endl;
      }
      return false;
    }
  }
  return true;
}

void SingleGl1PoolInputv2::ClearCurrentEvent()
{
  // called interactively, to get rid of the current event
  uint64_t currentbclk = *m_BclkStack.begin();
  //  std::cout << PHWHERE << "clearing bclk 0x" << std::hex << currentbclk << std::dec << std::endl;
  CleanupUsedPackets(currentbclk);
  // m_BclkStack.erase(currentbclk);
  return;
}

bool SingleGl1PoolInputv2::GetSomeMoreEvents()
{
  if (AllDone())
  {
    return false;
  }
  if (m_Gl1RawHitMap.empty())
  {
    return true;
  }

  uint64_t lowest_bclk = m_Gl1RawHitMap.begin()->first;
  lowest_bclk += m_BcoRange;
  uint64_t last_bclk = m_Gl1RawHitMap.rbegin()->first;
  if (Verbosity() > 1)
  {
    std::cout << PHWHERE << "first bclk 0x" << std::hex << lowest_bclk
              << " last bco: 0x" << last_bclk
              << std::dec << std::endl;
  }
  if (lowest_bclk >= last_bclk)
  {
    return true;
  }
  return false;
}

void SingleGl1PoolInputv2::CreateDSTNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  PHNodeIterator iterDst(dstNode);
  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", "GL1"));
  if (!detNode)
  {
    detNode = new PHCompositeNode("GL1");
    dstNode->addNode(detNode);
  }
  Gl1Packet *gl1hitcont = findNode::getClass<Gl1Packet>(detNode, "GL1RAWHIT");
  if (!gl1hitcont)
  {
    gl1hitcont = new Gl1Packetv2();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(gl1hitcont, "GL1RAWHIT", "PHObject");
    detNode->addNode(newNode);
  }
}
