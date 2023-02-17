#include "DumpTowerInfoContainer.h"

#include <phool/PHIODataNode.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<TowerInfoContainer>;

DumpTowerInfoContainer::DumpTowerInfoContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpTowerInfoContainer::process_Node(PHNode *myNode)
{
  TowerInfoContainer *towerinfocontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    towerinfocontainer = thisNode->getData();
  }
  if (towerinfocontainer)
  {
    *fout << "size: " << towerinfocontainer->size() << std::endl;
    for (auto hiter = towerinfocontainer->begin(); hiter !=towerinfocontainer->end(); ++hiter )
    {
      TowerInfo *rawtwr = hiter->second;
      *fout << "time: " << rawtwr->get_time() << std::endl;
      *fout << "energy: " << rawtwr->get_energy() << std::endl;
    }
  }
  return 0;
}
