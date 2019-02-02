#include "DumpRawTowerContainer.h"

#include <phool/PHIODataNode.h>

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>

#include <string>

using namespace std;

typedef PHIODataNode<RawTowerContainer> MyNode_t;

DumpRawTowerContainer::DumpRawTowerContainer(const string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpRawTowerContainer::process_Node(PHNode *myNode)
{
  RawTowerContainer *rawtowercontainer = NULL;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    rawtowercontainer = thisNode->getData();
  }
  if (rawtowercontainer)
  {
    RawTowerContainer::ConstIterator hiter;
    RawTowerContainer::ConstRange begin_end = rawtowercontainer->getTowers();
    *fout << "size: " << rawtowercontainer->size() << endl;
    for (hiter = begin_end.first; hiter != begin_end.second; ++hiter)
    {
      RawTower *rawtwr = hiter->second;
      *fout << "bineta: " << rawtwr->get_bineta() << endl;
      *fout << "binphi: " << rawtwr->get_binphi() << endl;
      *fout << "energy: " << rawtwr->get_energy() << endl;
      RawTower::CellConstRange cbegin_end = rawtwr->get_g4cells();
      for (RawTower::CellConstIterator iter = cbegin_end.first; iter != cbegin_end.second; ++iter)
      {
        *fout << "cell key: 0x" << hex << iter->first
              << dec << " edep: " << iter->second << endl;
      }
      RawTower::ShowerConstRange sbegin_end = rawtwr->get_g4showers();
      for (RawTower::ShowerConstIterator iter = sbegin_end.first; iter != sbegin_end.second; ++iter)
      {
        *fout << "shower id: " << iter->first
              << " edep: " << iter->second << endl;
      }
    }
  }
  return 0;
}
