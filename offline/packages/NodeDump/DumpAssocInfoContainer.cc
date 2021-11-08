#include "DumpAssocInfoContainer.h"

#include <phool/PHIODataNode.h>

#include <trackreco/AssocInfoContainer.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using namespace std;

typedef PHIODataNode<AssocInfoContainer> MyNode_t;

DumpAssocInfoContainer::DumpAssocInfoContainer(const string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpAssocInfoContainer::process_Node(PHNode *myNode)
{
  AssocInfoContainer *associnfocontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    associnfocontainer = thisNode->getData();
  }
  if (associnfocontainer)
  {
    associnfocontainer->identify(*fout);
  }
  return 0;
}
