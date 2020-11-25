#include "DumpRunHeader.h"

#include "DumpObject.h"

#include <ffaobjects/RunHeader.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using namespace std;

typedef PHIODataNode<RunHeader> MyNode_t;

DumpRunHeader::DumpRunHeader(const string &NodeName)
  : DumpObject(NodeName)
{
  WriteRunEvent(0);  // do not write info for each event
  node_written = 0;  // write runwise nodes only once
  return;
}

int DumpRunHeader::process_Node(PHNode *myNode)
{
  RunHeader *runheader = 0;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    runheader = thisNode->getData();
  }
  if (!node_written && runheader)
  {
    *fout << "RunHeader->isValid(): " << runheader->isValid() << endl;
    if (runheader->isValid())
    {
      *fout << "get_RunNumber(): " << runheader->get_RunNumber() << endl;
      node_written = 1;
    }
  }
  return 0;
}
