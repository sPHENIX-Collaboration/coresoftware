#include "SingleZdcTriggeredInput.h"

#include <frog/FROG.h>

#include <ffarawobjects/CaloPacketContainerv1.h>
#include <ffarawobjects/CaloPacketv1.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>
#include <Event/packet.h>

#include <TSystem.h>

#include <cstdint>   // for uint64_t
#include <iostream>  // for operator<<, basic_ostream, endl
#include <set>
#include <utility>  // for pair
#include <vector>

SingleZdcTriggeredInput::SingleZdcTriggeredInput(const std::string &name)
  : SingleTriggeredInput(name)
{
}

void SingleZdcTriggeredInput::CreateDSTNode(PHCompositeNode *my_topNode)
{
  std::array<std::string, 2> detname{"ZDC", "SEPD"};
  topNode(my_topNode);
  PHNodeIterator iter(topNode());
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode()->addNode(dstNode);
  }
  for (const auto &striter : detname)
  {
    PHNodeIterator iterDst(dstNode);
    PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", striter));
    if (!detNode)
    {
      detNode = new PHCompositeNode(striter);
      dstNode->addNode(detNode);
    }
    std::string OutNodeName = striter + std::string("Packets");
    CaloPacketContainer *packetcont = findNode::getClass<CaloPacketContainer>(detNode, OutNodeName);
    if (!packetcont)
    {
      packetcont = new CaloPacketContainerv1();
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(packetcont, OutNodeName, "PHObject");
      detNode->addNode(newNode);
    }
  }
}

void SingleZdcTriggeredInput::AddPacket(PHCompositeNode *topNode, OfflinePacket *newhit)
{
  std::string OutNodeName;
  if (std::clamp(newhit->getIdentifier(), 9000, 9999))  // sepd packet?
  {
    OutNodeName = "SEPD";
  }
  else
  {
    OutNodeName = "ZDC";
  }
  OutNodeName = OutNodeName + "Packets";
  CaloPacketContainer *packetcont = findNode::getClass<CaloPacketContainer>(topNode, OutNodeName);
  if (!packetcont)
  {
    std::cout << PHWHERE << " Could not locate " << OutNodeName << ", or type mismatch" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  CaloPacket *calopacket = dynamic_cast<CaloPacket *>(newhit);
  if (!calopacket)
  {
    std::cout << PHWHERE << " dynamic cast to CaloPacket failed for " << std::endl;
    newhit->identify();
    gSystem->Exit(1);
    exit(1);
  }
  packetcont->AddPacket(calopacket);
}
