#include "DumpTrkrClusterCrossingAssoc.h"

#include <phool/PHIODataNode.h>

#include <trackbase/TrkrClusterCrossingAssoc.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<TrkrClusterCrossingAssoc>;

DumpTrkrClusterCrossingAssoc::DumpTrkrClusterCrossingAssoc(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpTrkrClusterCrossingAssoc::process_Node(PHNode *myNode)
{
  TrkrClusterCrossingAssoc *trkrclustercrossingassoc = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    trkrclustercrossingassoc = thisNode->getData();
  }
  if (trkrclustercrossingassoc)
  {
    TrkrClusterCrossingAssoc::ConstRange begin_end = trkrclustercrossingassoc->getAll();
    *fout << "size " << trkrclustercrossingassoc->size() << std::endl;
    for (auto iter = begin_end.first; iter != begin_end.second; ++iter)
    {
      *fout << "cluster: " << std::hex << iter->first << std::dec << ", crossing: " << iter->second << std::endl;
    }
  }
  return 0;
}
