#include "DumpTowerBackground.h"

#include <jetbackground/TowerBackground.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<TowerBackground>;

DumpTowerBackground::DumpTowerBackground(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpTowerBackground::process_Node(PHNode *myNode)
{
  TowerBackground *twrbkg = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    twrbkg = thisNode->getData();
  }
  if (twrbkg)
  {
    *fout << "get_v2(): " << twrbkg->get_v2() << std::endl;
    *fout << "get_Psi2(): " << twrbkg->get_Psi2() << std::endl;
  }
  return 0;
}
