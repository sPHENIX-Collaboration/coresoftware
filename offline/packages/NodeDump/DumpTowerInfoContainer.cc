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
      *fout << "time_float: " << rawtwr->get_time_float() << std::endl;
      *fout << "chi2: " << rawtwr->get_chi2() << std::endl;
      *fout << "pedestal: " << rawtwr->get_pedestal() << std::endl;
      *fout << "isHot: " << rawtwr->get_isHot() << std::endl;
      *fout << "isBadTime: " << rawtwr->get_isBadTime() << std::endl;
      *fout << "isNotInstr: " << rawtwr->get_isNotInstr() << std::endl;
      *fout << "isGood: " << rawtwr->get_isGood() << std::endl;
      *fout << "status: " << rawtwr->get_status() << std::endl;
      *fout << "nsample: " << rawtwr->get_nsample() << std::endl;
      for (int j = 0; j < rawtwr->get_nsample(); j++)
      {
        *fout << "waveform_value[" << j << "]: " << rawtwr->get_waveform_value(j) << std::endl;
      }
    }
  }
  return 0;
}
