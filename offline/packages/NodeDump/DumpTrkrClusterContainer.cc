#include "DumpTrkrClusterContainer.h"

#include <phool/PHIODataNode.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

typedef PHIODataNode<TrkrClusterContainer> MyNode_t;

DumpTrkrClusterContainer::DumpTrkrClusterContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpTrkrClusterContainer::process_Node(PHNode *myNode)
{
  TrkrClusterContainer *trkrclustercontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    trkrclustercontainer = thisNode->getData();
  }
  if (trkrclustercontainer)
  {
    TrkrClusterContainer::ConstIterator hiter;
    *fout << "size: " << trkrclustercontainer->size() << std::endl;
    trkrclustercontainer->identify(*fout);
    TrkrClusterContainer::HitSetKeyList keylist = trkrclustercontainer->getHitSetKeys();
    for (auto iter = keylist.begin(); iter != keylist.end(); ++iter)
    {
      TrkrClusterContainer::ConstRange begin_end = trkrclustercontainer->getClusters(*iter);
      for (hiter = begin_end.first; hiter != begin_end.second; ++hiter)
      {
        TrkrCluster *trkrcluster = hiter->second;
        *fout << "getAdc: " << trkrcluster->getAdc() << std::endl;
        *fout << "getRPhiError: " << trkrcluster->getRPhiError() << std::endl;
        *fout << "getZError: " << trkrcluster->getZError() << std::endl;
      }
    }
  }
  return 0;
}
