#include "DumpTrkrClusterHitAssoc.h"

#include <phool/PHIODataNode.h>

#include <trackbase/TrkrClusterHitAssoc.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<TrkrClusterHitAssoc>;

DumpTrkrClusterHitAssoc::DumpTrkrClusterHitAssoc(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpTrkrClusterHitAssoc::process_Node(PHNode *myNode)
{
  TrkrClusterHitAssoc *trkrclusterhitassoc = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    trkrclusterhitassoc = thisNode->getData();
  }
  if (trkrclusterhitassoc)
  {
    trkrclusterhitassoc->identify(*fout);
  }
  return 0;
}
