#include "DumpTowerBackground.h"

#include <jetbackground/TowerBackground.h>

#include <phool/PHIODataNode.h>

#include <string>


using namespace std;

typedef PHIODataNode<TowerBackground> MyNode_t;

DumpTowerBackground::DumpTowerBackground(const string &NodeName): DumpObject(NodeName)
{
  return ;
}

int DumpTowerBackground::process_Node(PHNode *myNode)
{
  TowerBackground *twrbkg = NULL;
  MyNode_t *thisNode = static_cast <MyNode_t *> (myNode);
  if (thisNode)
  {
    twrbkg = thisNode->getData();
  }
  if (twrbkg)
  {
    *fout << "get_v2(): "  << twrbkg->get_v2() << endl;
    *fout << "get_Psi2(): "   << twrbkg->get_Psi2() << endl;
  }
  return 0;
}
