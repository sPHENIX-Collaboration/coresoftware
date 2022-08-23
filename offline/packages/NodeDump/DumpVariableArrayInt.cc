#include "DumpVariableArrayInt.h"

#include <phool/PHIODataNode.h>
#include <vararray/VariableArrayInt.h>

#include <string>

#include <cstddef>
#include <map>
#include <ostream>
#include <utility>

typedef PHIODataNode<VariableArrayInt> MyNode_t;

DumpVariableArrayInt::DumpVariableArrayInt(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpVariableArrayInt::process_Node(PHNode *myNode)
{
  VariableArrayInt *variablearray = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    variablearray = thisNode->getData();
  }
  if (variablearray)
  {
    *fout << "Id(): " << variablearray->Id() << std::endl;
    *fout << "get_array_size(): " << variablearray->get_array_size() << std::endl;
    const int *sval = variablearray->get_array();
    for (unsigned int i = 0; i < variablearray->get_array_size(); i++)
    {
      *fout << "val[" << i << "]: " << sval[i] << std::endl;
    }
  }
  return 0;
}
