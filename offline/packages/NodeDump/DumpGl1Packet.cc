#include "DumpGl1Packet.h"

#include <ffarawobjects/Gl1Packet.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<Gl1Packet>;

DumpGl1Packet::DumpGl1Packet(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpGl1Packet::process_Node(PHNode *myNode)
{
  Gl1Packet *gl1packet = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    gl1packet = thisNode->getData();
  }
  if (gl1packet)
  {
    *fout << "packet_nr: " << gl1packet->iValue(0) << std::endl;
    *fout << "BCO: 0x" << gl1packet->lValue(0, "BCO") << std::endl;
    *fout << "TriggerInput: " << gl1packet->lValue(0, "TriggerInput") << std::endl;
    *fout << "TriggerVector: " << gl1packet->lValue(0, "TriggerVector") << std::endl;
    *fout << "BunchNumber: " << gl1packet->lValue(0, "BunchNumber") << std::endl;
    for (int i = 0; i < 64; i++)
    {
      *fout << "lValue(" << i << ",0): " << gl1packet->lValue(i, 0) << std::endl;
      *fout << "lValue(" << i << ",1): " << gl1packet->lValue(i, 1) << std::endl;
      *fout << "lValue(" << i << ",2): " << gl1packet->lValue(i, 2) << std::endl;
    }
  }
  return 0;
}
