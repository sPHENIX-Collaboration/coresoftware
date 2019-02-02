#include "DumpSyncObject.h"

#include <ffaobjects/SyncObject.h>

#include <phool/PHIODataNode.h>

#include <string>

using namespace std;

typedef PHIODataNode<SyncObject> MyNode_t;

DumpSyncObject::DumpSyncObject(const string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpSyncObject::process_Node(PHNode *myNode)
{
  SyncObject *syncobject = NULL;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    syncobject = thisNode->getData();
  }
  if (syncobject)
  {
    *fout << "SyncObject->isValid(): " << syncobject->isValid() << endl;
    if (syncobject->isValid())
    {
      *fout << "EventCounter(): " << syncobject->EventCounter() << endl;
      *fout << "EventNumber(): " << syncobject->EventNumber() << endl;
      *fout << "RunNumber(): " << syncobject->RunNumber() << endl;
      *fout << "SegmentNumber(): " << syncobject->SegmentNumber() << endl;
    }
  }
  return 0;
}
