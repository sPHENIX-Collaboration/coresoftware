#include "DumpRunHeader.h"

#include "DumpObject.h"

#include <ffaobjects/RunHeader.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<RunHeader>;

DumpRunHeader::DumpRunHeader(const std::string &NodeName)
  : DumpObject(NodeName)
{
  WriteRunEvent(0);  // do not write info for each event
  return;
}

int DumpRunHeader::process_Node(PHNode *myNode)
{
  RunHeader *runheader = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    runheader = thisNode->getData();
  }
  if (!node_written && runheader)
  {
    *fout << "RunHeader->isValid(): " << runheader->isValid() << std::endl;
    if (runheader->isValid())
    {
      *fout << "get_RunNumber(): " << runheader->get_RunNumber() << std::endl;
      node_written = 1;
    }
  }
  return 0;
}
