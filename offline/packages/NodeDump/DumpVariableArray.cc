#include "DumpVariableArray.h"

#include <vararray/VariableArray.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<VariableArray>;

DumpVariableArray::DumpVariableArray(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpVariableArray::process_Node(PHNode *myNode)
{
  VariableArray *variablearray = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    variablearray = thisNode->getData();
  }
  if (variablearray)
  {
    *fout << "Id(): " << variablearray->Id() << std::endl;
    *fout << "get_array_size(): " << variablearray->get_array_size() << std::endl;
    const short *sval = variablearray->get_array();
    for (unsigned int i = 0; i < variablearray->get_array_size(); i++)
    {
      *fout << "val[" << i << "]: " << sval[i] << std::endl;
    }
  }
  return 0;
}
