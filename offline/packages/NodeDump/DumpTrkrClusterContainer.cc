#include "DumpTrkrClusterContainer.h"

#include <phool/PHIODataNode.h>

#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <cmath>
#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<TrkrClusterContainer>;

DumpTrkrClusterContainer::DumpTrkrClusterContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpTrkrClusterContainer::process_Node(PHNode *myNode)
{
  TrkrClusterContainer *trkrclustercontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);  // NOLINT(cppcoreguidelines-pro-type-static-cast-downcast)
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
    for (unsigned int &iter : keylist)
    {
      TrkrClusterContainer::ConstRange begin_end = trkrclustercontainer->getClusters(iter);
      for (hiter = begin_end.first; hiter != begin_end.second; ++hiter)
      {
        const TrkrDefs::cluskey cluster_key = hiter->first;
        TrkrCluster *trkrcluster = hiter->second;
        const auto cluster_errors = ClusterErrorPara::get_clusterv5_modified_error(trkrcluster, 0.0, cluster_key);
        *fout << "getAdc: " << trkrcluster->getAdc() << std::endl;
        *fout << "getRPhiError: " << trkrcluster->getRPhiError() << std::endl;
        *fout << "getZError: " << trkrcluster->getZError() << std::endl;
        *fout << "ClusterErrorPara RPhiError: " << std::sqrt(cluster_errors.first) << std::endl;
        *fout << "ClusterErrorPara ZError: " << std::sqrt(cluster_errors.second) << std::endl;
      }
    }
  }
  return 0;
}
