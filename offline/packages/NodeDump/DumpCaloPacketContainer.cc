#include "DumpCaloPacketContainer.h"

#include <ffarawobjects/CaloPacket.h>
#include <ffarawobjects/CaloPacketContainer.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<CaloPacketContainer>;

DumpCaloPacketContainer::DumpCaloPacketContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpCaloPacketContainer::process_Node(PHNode *myNode)
{
  CaloPacketContainer *calocont = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    calocont = thisNode->getData();
  }
  if (calocont)
  {
    for (unsigned int i = 0; i < calocont->get_npackets(); i++)
    {
      CaloPacket *calopacket = calocont->getPacket(i);
      *fout << "packet_nr: " << calopacket->iValue(0) << std::endl;
      *fout << "EventNr: " << calopacket->iValue(0, "EVTNR") << std::endl;
      *fout << "Clock: 0x" << std::hex << calopacket->iValue(0, "CLOCK") << std::dec << std::endl;
      *fout << "Modules: " << calopacket->iValue(0, "NRMODULES") << std::endl;
      *fout << "Channels: " << calopacket->iValue(0, "CHANNELS") << std::endl;
      *fout << "Samples: " << calopacket->iValue(0, "SAMPLES") << std::endl;
      *fout << "DetId: 0x" << std::hex << calopacket->iValue(0, "DETID") << std::dec << std::endl;
      *fout << "Samples: " << std::hex << calopacket->iValue(0, "MODULEADDRESS") << std::dec << std::endl;
      *fout << "EVENCHECKSUMOK: " << calopacket->iValue(0, "EVENCHECKSUMOK") << std::endl;
      *fout << "ODDCHECKSUMOK: " << calopacket->iValue(0, "ODDCHECKSUMOK") << std::endl;
      for (int k = 0; k < calopacket->iValue(0, "NRMODULES"); k++)
      {
        *fout << "FEM Slot: " << calopacket->iValue(k, "FEMSLOT") << std::endl;
        *fout << "FEM Evt nr: " << calopacket->iValue(k, "FEMEVTNR") << std::endl;
        *fout << "FEM Clock: " << calopacket->iValue(k, "FEMCLOCK") << std::endl;
        *fout << "FEM Checksum LSB: 0x" << std::hex << calopacket->iValue(k, "CHECKSUMLSB") << std::dec << std::endl;
        *fout << "FEM Checksum MSB: 0x" << std::hex << calopacket->iValue(k, "CHECKSUMMSB") << std::dec << std::endl;
      }
      for (int k = 0; k < calopacket->iValue(0, "CHANNELS"); k++)
      {
        *fout << "SUPPRESSED: " << calopacket->iValue(k, "SUPPRESSED") << std::endl;
        *fout << "PRE: 0x" << std::hex << calopacket->iValue(k, "PRE") << std::dec << std::endl;
        *fout << "POST: 0x" << std::hex << calopacket->iValue(k, "POST") << std::dec << std::endl;
        for (int j = 0; j < calopacket->iValue(0, "SAMPLES"); j++)
        {
          *fout << "iValue(" << j << ", " << k << "): " << calopacket->iValue(j, k) << std::endl;
        }
      }
    }
  }
  return 0;
}
