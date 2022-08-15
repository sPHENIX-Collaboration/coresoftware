#include "DumpCdbUrlSave.h"

#include <phool/PHIODataNode.h>

#include <ffaobjects/CdbUrlSave.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

typedef PHIODataNode<CdbUrlSave> MyNode_t;

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
