#include "DumpTpcSeedTrackMap.h"

#include <phool/PHIODataNode.h>

#include <trackbase/TpcSeedTrackMap.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<TpcSeedTrackMap>;

DumpTpcSeedTrackMap::DumpTpcSeedTrackMap(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpTpcSeedTrackMap::process_Node(PHNode *myNode)
{
  TpcSeedTrackMap *tpcseedtrackmap = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    tpcseedtrackmap = thisNode->getData();
  }
  if (tpcseedtrackmap)
  {
    TpcSeedTrackMap::ConstRange begin_end = tpcseedtrackmap->getAll();
    *fout << "size " << tpcseedtrackmap->size() << std::endl;
    for (auto iter = begin_end.first; iter != begin_end.second; ++iter)
    {
      *fout << "original id: " << iter->first << ", duplicate id: " << iter->second << std::endl;
    }
  }
  return 0;
}
