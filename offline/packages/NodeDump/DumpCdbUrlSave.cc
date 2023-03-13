#include "DumpCdbUrlSave.h"

#include <phool/PHIODataNode.h>

#include <ffaobjects/CdbUrlSave.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<CdbUrlSave>;

DumpCdbUrlSave::DumpCdbUrlSave(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpCdbUrlSave::process_Node(PHNode *myNode)
{
  CdbUrlSave *cdburlsave = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    cdburlsave = thisNode->getData();
  }
  if (cdburlsave)
  {
    cdburlsave->identify(*fout);
  }
  return 0;
}
