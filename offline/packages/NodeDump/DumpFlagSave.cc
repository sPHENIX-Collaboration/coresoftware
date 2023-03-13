#include "DumpFlagSave.h"

#include <phool/PHIODataNode.h>

#include <ffaobjects/FlagSave.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<FlagSave>;

DumpFlagSave::DumpFlagSave(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpFlagSave::process_Node(PHNode *myNode)
{
  FlagSave *flagsave = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    flagsave = thisNode->getData();
  }
  if (flagsave)
  {
    flagsave->identify(*fout);
  }
  return 0;
}
