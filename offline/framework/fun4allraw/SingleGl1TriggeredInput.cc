#include "SingleGl1TriggeredInput.h"

#include "Fun4AllPrdfInputTriggerManager.h"
#include "InputManagerType.h"

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

#include <TSystem.h>

#include <cstdint>   // for uint64_t
#include <iostream>  // for operator<<, basic_ostream<...
#include <iterator>  // for reverse_iterator
#include <limits>    // for numeric_limits
#include <memory>
#include <set>
#include <utility>  // for pair

SingleGl1TriggeredInput::SingleGl1TriggeredInput(const std::string &name)
  : SingleTriggeredInput(name)
{
}

SingleGl1TriggeredInput::~SingleGl1TriggeredInput()
{
}

void SingleGl1TriggeredInput::FillPool(const unsigned int keep)
{
  if (AllDone())  // no more files and all events read
  {
    return;
  }
//  OfflinePacket *gl1hitcont = findNode::getClass<OfflinePacket>(m_topNode, "GL1Packet");
  if (!FilesDone())
  {
  FillEventVector();
  }
  if (keep > 100000000)
  {
    std::cout << "huh?" << std::endl;
  }
  if (m_EventDeque.empty())
  {
    std::cout << Name() << ":all events done" << std::endl;
    AllDone(1);
    return;
  }
  Event *evt = m_EventDeque.front();
  m_EventDeque.pop_front();
  RunNumber(evt->getRunNumber());
    int EventSequence = evt->getEvtSequence();
    if (keep > 100000000)
    {
    std::cout << "run: " << RunNumber() << ", evt: " << EventSequence << std::endl;
    }
    Packet *packet = evt->getPacket(14001);
    if (packet)
    {
    Gl1Packet *newhit = new Gl1Packetv2();
    uint64_t gtm_bco = packet->lValue(0, "BCO");
//    std::cout << "saving bco 0x" << std::hex << gtm_bco << std::dec << std::endl;
    unsigned int packetnumber = packet->iValue(0);
    newhit->setBCO(packet->lValue(0, "BCO"));
    newhit->setHitFormat(packet->getHitFormat());
    newhit->setIdentifier(packet->getIdentifier());
    newhit->setEvtSequence(EventSequence);
    newhit->setPacketNumber(packetnumber);

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
    delete packet;
       delete newhit;
    }
       delete evt;
    
}


void SingleGl1TriggeredInput::Print(const std::string &what) const
{
  std::cout << "what: " << what << std::endl;
}

void SingleGl1TriggeredInput::CreateDSTNode(PHCompositeNode *topNode)
{
  m_topNode = topNode;
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
  OfflinePacket *gl1hitcont = findNode::getClass<OfflinePacket>(detNode, "GL1Packet");
  if (!gl1hitcont)
  {
    gl1hitcont = new Gl1Packetv2();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(gl1hitcont, "GL1Packet", "PHObject");
    detNode->addNode(newNode);
  }
}

uint64_t SingleGl1TriggeredInput::GetClock(Event *evt)
{
  Packet *packet = evt->getPacket(14001);
uint64_t clock = packet->lValue(0, "BCO");
delete packet;
return clock;
}
