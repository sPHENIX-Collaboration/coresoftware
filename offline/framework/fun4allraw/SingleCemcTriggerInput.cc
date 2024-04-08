#include "SingleCemcTriggerInput.h"

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

static const int NCEMCPACKETS = 17;

SingleCemcTriggerInput::SingleCemcTriggerInput(const std::string &name)
  : SingleTriggerInput(name)
{
  SubsystemEnum(InputManagerType::CEMC);
  plist = new Packet *[NCEMCPACKETS];  // 16 packets for the cemc in each file
}

SingleCemcTriggerInput::~SingleCemcTriggerInput()
{
  CleanupUsedPackets(std::numeric_limits<int>::max());
  delete[] plist;
}

void SingleCemcTriggerInput::FillPool(const unsigned int /*nbclks*/)
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
    int npackets = evt->getPacketList(plist, NCEMCPACKETS); // just in case we have more packets, they will not vanish silently
    if (npackets >= NCEMCPACKETS)
    {
      std::cout << PHWHERE << " Packets array size " << NCEMCPACKETS
		<< " too small for " << Name()
<< ", increase NCEMCPACKETS and rebuild" << std::endl;
      exit(1);
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
      int nr_channels = plist[i]->iValue(0, "CHANNELS");
      int nr_samples = plist[i]->iValue(0, "SAMPLES");
      if (nr_modules > newhit->getMaxNumModules() )
      {
	std::cout << PHWHERE << " too many modules " <<  nr_modules << ", max is " 
                  <<  newhit->getMaxNumModules() << ", need to adjust arrays" << std::endl;
	gSystem->Exit(1);
      }
      if (nr_channels > newhit->getMaxNumChannels())
      {
	std::cout << PHWHERE << " too many channels " <<  nr_channels << ", max is " 
                  <<  newhit->getMaxNumChannels() << ", need to adjust arrays" << std::endl;
	gSystem->Exit(1);
      }
      if (nr_samples > newhit->getMaxNumSamples())
      {
	std::cout << PHWHERE << " too many samples " <<  nr_samples << ", max is " 
                  <<  newhit->getMaxNumSamples() << ", need to adjust arrays" << std::endl;
	gSystem->Exit(1);
      }

      uint64_t gtm_bco = plist[i]->iValue(0, "CLOCK");
      newhit->setNrModules(nr_modules);
      newhit->setNrSamples(nr_samples);
      newhit->setNrChannels(nr_channels);
      newhit->setBCO(gtm_bco);
      newhit->setPacketEvtSequence(plist[i]->iValue(0, "EVTNR"));
      newhit->setIdentifier(plist[i]->getIdentifier());
      newhit->setHitFormat(plist[i]->getHitFormat());
      newhit->setEvtSequence(EventSequence);
      newhit->setEvenChecksum(plist[i]->iValue(0, "EVENCHECKSUM"));
      newhit->setCalcEvenChecksum(plist[i]->iValue(0, "CALCEVENCHECKSUM"));
      newhit->setOddChecksum(plist[i]->iValue(0, "ODDCHECKSUM"));
      newhit->setCalcOddChecksum(plist[i]->iValue(0, "CALCODDCHECKSUM"));
      newhit->setModuleAddress(plist[i]->iValue(0,"MODULEADDRESS"));
      newhit->setDetId(plist[i]->iValue(0,"DETID"));
      for (int ifem = 0; ifem < nr_modules; ifem++)
      {
        newhit->setFemClock(ifem, plist[i]->iValue(ifem, "FEMCLOCK"));
        newhit->setFemEvtSequence(ifem, plist[i]->iValue(ifem, "FEMEVTNR"));
        newhit->setFemSlot(ifem, plist[i]->iValue(ifem, "FEMSLOT"));
      }
      for (int ipmt = 0; ipmt < nr_channels; ipmt++)
      {
        for (int isamp = 0; isamp < nr_samples; isamp++)
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
        TriggerInputManager()->AddCemcPacket(EventSequence, newhit);
      }
      m_CemcPacketMap[EventSequence].push_back(newhit);
      m_EventStack.insert(EventSequence);
      if (ddump_enabled())
      {
	ddumppacket(plist[i]);
      }
      delete plist[i];
    }
  }
}

void SingleCemcTriggerInput::Print(const std::string &what) const
{
  if (what == "ALL" || what == "STORAGE")
  {
    for (const auto &bcliter : m_CemcPacketMap)
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

void SingleCemcTriggerInput::CleanupUsedPackets(const int eventno)
{
  std::vector<int> toclearevents;
  for (const auto &iter : m_CemcPacketMap)
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
    m_CemcPacketMap.erase(iter);
  }
}

void SingleCemcTriggerInput::ClearCurrentEvent()
{
  // called interactively, to get rid of the current event
  int currentevent = *m_EventStack.begin();
  //  std::cout << PHWHERE << "clearing bclk 0x" << std::hex << currentbclk << std::dec << std::endl;
  CleanupUsedPackets(currentevent);
  return;
}

bool SingleCemcTriggerInput::GetSomeMoreEvents()
{
  if (AllDone())
  {
    return false;
  }
  if (m_CemcPacketMap.empty())
  {
    return true;
  }

  int first_event = m_CemcPacketMap.begin()->first;
  int last_event = m_CemcPacketMap.rbegin()->first;
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

void SingleCemcTriggerInput::CreateDSTNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  PHNodeIterator iterDst(dstNode);
  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", "CEMC"));
  if (!detNode)
  {
    detNode = new PHCompositeNode("CEMC");
    dstNode->addNode(detNode);
  }
  CaloPacketContainer *cemcpacketcont = findNode::getClass<CaloPacketContainer>(detNode, "CEMCPackets");
  if (!cemcpacketcont)
  {
    cemcpacketcont = new CaloPacketContainerv1();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(cemcpacketcont, "CEMCPackets", "PHObject");
    detNode->addNode(newNode);
  }
}
