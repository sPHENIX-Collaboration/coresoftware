#include "SingleHcalTriggerInput.h"

#include "Fun4AllPrdfInputTriggerManager.h"
#include "InputManagerType.h"

#include <ffarawobjects/CaloPacketv1.h>
#include <ffarawobjects/CaloPacketContainerv1.h>

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

#include <TSystem.h>

#include <cstdint>   // for uint64_t
#include <iostream>  // for operator<<, basic_ostream<...
#include <iterator>  // for reverse_iterator
#include <limits>    // for numeric_limits
#include <memory>
#include <set>
#include <utility>  // for pair

SingleHcalTriggerInput::SingleHcalTriggerInput(const std::string &name)
  : SingleTriggerInput(name)
{
  SubsystemEnum(InputManagerType::HCAL);
  plist = new Packet *[2];  // two packets for the hcal in each file
}

SingleHcalTriggerInput::~SingleHcalTriggerInput()
{
  CleanupUsedPackets(std::numeric_limits<int>::max());
  delete[] plist;
}

void SingleHcalTriggerInput::FillPool(const unsigned int /*nbclks*/)
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
    RunNumber(evt->getRunNumber());
    if (GetVerbosity() > 1)
    {
      evt->identify();
    }
    if (evt->getEvtType() != DATAEVENT)
    {
      m_NumSpecialEvents++;
      continue;
    }
    int EventSequence = evt->getEvtSequence();
    int npackets = evt->getPacketList(plist, 2);
    if (npackets > 2)
    {
      std::cout << PHWHERE << " Number of packets in array (2) too small for "
                << npackets << std::endl;
      exit(1);
    }

    for (int i = 0; i < npackets; i++)
    {
      if (Verbosity() > 2)
      {
        plist[i]->identify();
      }

      // by default use previous bco clock for gtm bco
      CaloPacket *newhit = new CaloPacketv1();
      int nr_modules = plist[i]->iValue(0,"NRMODULES");
      if (nr_modules > 3)
      {
	std::cout << PHWHERE << " too many modules, need to adjust arrays" << std::endl;
	gSystem->Exit(1);
      }
      uint64_t gtm_bco = plist[i]->iValue(0, "CLOCK");
      newhit->setBCO(plist[i]->iValue(0, "CLOCK"));
      newhit->setPacketEvtSequence(plist[i]->iValue(0, "EVTNR"));
      newhit->setIdentifier(plist[i]->getIdentifier());
      newhit->setEvtSequence(EventSequence);
      for (int ifem = 0; ifem < nr_modules; ifem++)
      {
        newhit->setFemClock(ifem, plist[i]->iValue(ifem, "FEMCLOCK"));
      }
      for (int ipmt = 0; ipmt < 192; ipmt++)
      {
        for (int isamp = 0; isamp < 31; isamp++)
        {
          newhit->setSample(ipmt, isamp, plist[i]->iValue(isamp, ipmt));
        }
      }
      if (Verbosity() > 2)
      {
        std::cout << PHWHERE << "evtno: " << EventSequence
                  << ", bco: 0x" << std::hex << gtm_bco << std::dec
                  << std::endl;
      }
      if (TriggerInputManager())
      {
        TriggerInputManager()->AddHcalPacket(EventSequence, newhit);
      }
      m_HcalPacketMap[EventSequence].push_back(newhit);
      m_EventStack.insert(EventSequence);

      delete plist[i];
    }
  }
}

void SingleHcalTriggerInput::Print(const std::string &what) const
{
  if (what == "ALL" || what == "STORAGE")
  {
    for (const auto &bcliter : m_HcalPacketMap)
    {
      std::cout << PHWHERE << "Event: " << bcliter.first << std::endl;
    }
  }
  if (what == "ALL" || what == "STACK")
  {
    for (auto iter : m_EventStack)
    {
      std::cout << PHWHERE << "stacked event: " << iter << std::endl;
    }
  }
}

void SingleHcalTriggerInput::CleanupUsedPackets(const int eventno)
{
  std::vector<int> toclearevents;
  for (const auto &iter : m_HcalPacketMap)
  {
    if (iter.first <= eventno)
    {
      for (auto pktiter : iter.second)
      {
        delete pktiter;
      }
      toclearevents.push_back(iter.first);
    }
    else
    {
      break;
    }
  }

  for (auto iter : toclearevents)
  {
    m_EventStack.erase(iter);
    m_HcalPacketMap.erase(iter);
  }
}

void SingleHcalTriggerInput::ClearCurrentEvent()
{
  // called interactively, to get rid of the current event
  int currentevent = *m_EventStack.begin();
  //  std::cout << PHWHERE << "clearing bclk 0x" << std::hex << currentbclk << std::dec << std::endl;
  CleanupUsedPackets(currentevent);
  return;
}

bool SingleHcalTriggerInput::GetSomeMoreEvents()
{
  if (AllDone())
  {
    return false;
  }
  if (m_HcalPacketMap.empty())
  {
    return true;
  }

  int first_event = m_HcalPacketMap.begin()->first;
  int last_event = m_HcalPacketMap.rbegin()->first;
  if (Verbosity() > 1)
  {
    std::cout << PHWHERE << "first event: " << first_event
              << " last event: " << last_event
              << std::endl;
  }
  if (first_event >= last_event)
  {
    return true;
  }
  return false;
}

void SingleHcalTriggerInput::CreateDSTNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  PHNodeIterator iterDst(dstNode);
  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", "HCAL"));
  if (!detNode)
  {
    detNode = new PHCompositeNode("HCAL");
    dstNode->addNode(detNode);
  }
  CaloPacketContainer *hcalpacketcont = findNode::getClass<CaloPacketContainer>(detNode, "HCALPackets");
  if (!hcalpacketcont)
  {
    hcalpacketcont = new CaloPacketContainerv1();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(hcalpacketcont, "HCALPackets", "PHObject");
    detNode->addNode(newNode);
  }
}
