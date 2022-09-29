#include "DumpTrkrHitTruthClusters.h"

#include <trackbase/TrkrHitTruthClusters.h>

#include <phool/PHIODataNode.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<TrkrHitTruthClusters>;

DumpTrkrHitTruthClusters::DumpTrkrHitTruthClusters(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpTrkrHitTruthClusters::process_Node(PHNode *myNode)
{
  TrkrHitTruthClusters *trkrhittruthclusters = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    trkrhittruthclusters = thisNode->getData();
  }
  if (trkrhittruthclusters)
  {
    trkrhittruthclusters->identify(*fout);
  }
  return 0;
}
