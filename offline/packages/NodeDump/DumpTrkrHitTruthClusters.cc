#include "DumpTrkrHitTruthClusters.h"

#include <trackbase/TrkrHitTruthClusters.h>

#include <phool/PHIODataNode.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

typedef PHIODataNode<TrkrHitTruthClusters> MyNode_t;

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
    trkrhittruthclusters->print_clusters(*fout);
  }
  return 0;
}
