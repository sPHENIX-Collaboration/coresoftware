#include "DumpInttDeadMap.h"

#include <phool/PHIODataNode.h>

#include <g4intt/InttDeadMap.h>

#include <map>
#include <ostream>
#include <set>
#include <string>
#include <utility>

typedef PHIODataNode<InttDeadMap> MyNode_t;

DumpInttDeadMap::DumpInttDeadMap(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpInttDeadMap::process_Node(PHNode *myNode)
{
  InttDeadMap *inttdeadmap = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    inttdeadmap = thisNode->getData();
  }
  if (inttdeadmap)
  {
    *fout << "size " << inttdeadmap->size() << std::endl;
    const InttDeadMap::Map thismap = inttdeadmap->getDeadChannels();
    for (auto iter = thismap.begin(); iter != thismap.end(); ++iter)
    {
      *fout << "dead channel: " << std::hex << *iter << std::dec << std::endl;
    }
  }
  return 0;
}
