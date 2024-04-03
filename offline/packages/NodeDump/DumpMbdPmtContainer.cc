#include "DumpMbdPmtContainer.h"

#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<MbdPmtContainer>;

DumpMbdPmtContainer::DumpMbdPmtContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpMbdPmtContainer::process_Node(PHNode *myNode)
{
  MbdPmtContainer *mbdpmts = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    mbdpmts = thisNode->getData();
  }
  if (mbdpmts && mbdpmts->isValid())
  {
    *fout << "MbdPmtContainer->get_npmt: " << mbdpmts->get_npmt() << std::endl;
    for (int j = 0; j < mbdpmts->get_npmt(); j++)
    {
      *fout << "MbdPmtContainer->get_pmt(" << j << ")->get_q(): " << mbdpmts->get_pmt(j)->get_q() << std::endl;
      *fout << "MbdPmtContainer->get_pmt(" << j << ")->get_time(): " << mbdpmts->get_pmt(j)->get_time() << std::endl;
    }
  }
  return 0;
}
