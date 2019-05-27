#include "DumpRawClusterContainer.h"

#include <phool/PHIODataNode.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using namespace std;

typedef PHIODataNode<RawClusterContainer> MyNode_t;

DumpRawClusterContainer::DumpRawClusterContainer(const string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpRawClusterContainer::process_Node(PHNode *myNode)
{
  RawClusterContainer *rawclustercontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    rawclustercontainer = thisNode->getData();
  }
  if (rawclustercontainer)
  {
    RawClusterContainer::ConstIterator hiter;
    RawClusterContainer::ConstRange begin_end = rawclustercontainer->getClusters();
    *fout << "size: " << rawclustercontainer->size() << endl;
    for (hiter = begin_end.first; hiter != begin_end.second; ++hiter)
    {
      *fout << "NTowers: " << hiter->second->getNTowers() << endl;
      *fout << "z: " << hiter->second->get_z() << endl;
      *fout << "phi: " << hiter->second->get_phi() << endl;
      *fout << "energy: " << hiter->second->get_energy() << endl;
    }
  }
  return 0;
}
