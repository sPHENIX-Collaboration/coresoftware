#include "DumpTowerInfoContainer.h"

#include <phool/PHIODataNode.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <ostream>
#include <string>

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
    unsigned int nchannels = towerinfocontainer->size();
    *fout << "size: " << towerinfocontainer->size() << std::endl;
    for (unsigned int channel = 0; channel < nchannels; channel++)
    {
      TowerInfo *rawtwr = towerinfocontainer->get_tower_at_channel(channel);
      *fout << "time: " << rawtwr->get_time() << std::endl;
      *fout << "energy: " << rawtwr->get_energy() << std::endl;
    }
  }
  return 0;
}
