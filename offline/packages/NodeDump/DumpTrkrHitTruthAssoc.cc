#include "DumpTrkrHitTruthAssoc.h"

#include <phool/PHIODataNode.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<TrkrClusterContainer>;

DumpTrkrHitTruthAssoc::DumpTrkrHitTruthAssoc(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpTrkrHitTruthAssoc::process_Node(PHNode *myNode)
{
  TrkrClusterContainer *trkrclustercontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    trkrclustercontainer = thisNode->getData();
  }
  if (trkrclustercontainer)
  {
    trkrclustercontainer->identify(*fout);
  }
  return 0;
}
