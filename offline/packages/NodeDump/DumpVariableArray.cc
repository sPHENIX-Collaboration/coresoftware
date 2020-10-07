#include "DumpVariableArray.h"

#include <vararray/VariableArray.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using namespace std;

typedef PHIODataNode<VariableArray> MyNode_t;

DumpVariableArray::DumpVariableArray(const string &NodeName)
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
    *fout << "Id(): " << variablearray->Id() << endl;
    *fout << "get_array_size(): " << variablearray->get_array_size() << endl;
    const short *sval = variablearray->get_array();
    for (unsigned int i = 0; i < variablearray->get_array_size(); i++)
    {
      *fout << "val[" << i << "]: " << sval[i] << endl;
    }
  }
  return 0;
}
