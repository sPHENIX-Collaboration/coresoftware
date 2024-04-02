#include "DumpMvtxRawEvtHeader.h"

#include <ffarawobjects/MvtxRawEvtHeader.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<MvtxRawEvtHeader>;

DumpMvtxRawEvtHeader::DumpMvtxRawEvtHeader(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpMvtxRawEvtHeader::process_Node(PHNode *myNode)
{
  MvtxRawEvtHeader *mvtxRawEvtHeader = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    mvtxRawEvtHeader = thisNode->getData();
  }
  if (mvtxRawEvtHeader)
  {
    for (auto iter : mvtxRawEvtHeader->getMvtxFeeIdSet())
    {
      *fout << "FeeIdSet: " << iter << std::endl;
    }
    for (auto iter : mvtxRawEvtHeader->getMvtxLvL1BCO())
    {
      *fout << "LvL1BCO: " << iter << std::endl;
    }
  }
  return 0;
}
