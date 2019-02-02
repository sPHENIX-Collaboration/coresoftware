#include "DumpVariableArrayInt.h"

#include <phool/PHIODataNode.h>
#include <vararray/VariableArrayInt.h>

#include <string>

using namespace std;

typedef PHIODataNode<VariableArrayInt> MyNode_t;

DumpVariableArrayInt::DumpVariableArrayInt(const string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpVariableArrayInt::process_Node(PHNode *myNode)
{
  VariableArrayInt *variablearray = NULL;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    variablearray = thisNode->getData();
  }
  if (variablearray)
  {
    *fout << "Id(): " << variablearray->Id() << endl;
    *fout << "get_array_size(): " << variablearray->get_array_size() << endl;
    const int *sval = variablearray->get_array();
    for (unsigned int i = 0; i < variablearray->get_array_size(); i++)
    {
      *fout << "val[" << i << "]: " << sval[i] << endl;
    }
  }
  return 0;
}
