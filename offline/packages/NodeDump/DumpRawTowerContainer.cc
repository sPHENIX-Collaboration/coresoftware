#include "DumpRawTowerContainer.h"

#include <phool/PHIODataNode.h>

#include <g4cemc/RawTowerContainer.h>
#include <g4cemc/RawTower.h>

#include <string>

using namespace std;

typedef PHIODataNode<RawTowerContainer> MyNode_t;

DumpRawTowerContainer::DumpRawTowerContainer(const string &NodeName): DumpObject(NodeName)
{
  return ;
}

int DumpRawTowerContainer::process_Node(PHNode *myNode)
{
  RawTowerContainer *rawtowercontainer = NULL;
  MyNode_t *thisNode = static_cast <MyNode_t *> (myNode);
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
          *fout << "bineta: " << hiter->second->get_bineta() << endl;
          *fout << "binphi: " << hiter->second->get_binphi() << endl;
          *fout << "energy: " << hiter->second->get_energy() << endl;
        }
    }
  return 0;
}

