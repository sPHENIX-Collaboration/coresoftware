#include "DumpSyncObject.h"

#include <ffaobjects/SyncObject.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<SyncObject>;

DumpSyncObject::DumpSyncObject(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpSyncObject::process_Node(PHNode *myNode)
{
  SyncObject *syncobject = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    syncobject = thisNode->getData();
  }
  if (syncobject)
  {
    *fout << "SyncObject->isValid(): " << syncobject->isValid() << std::endl;
    if (syncobject->isValid())
    {
      *fout << "EventCounter(): " << syncobject->EventCounter() << std::endl;
      *fout << "EventNumber(): " << syncobject->EventNumber() << std::endl;
      *fout << "RunNumber(): " << syncobject->RunNumber() << std::endl;
      *fout << "SegmentNumber(): " << syncobject->SegmentNumber() << std::endl;
    }
  }
  return 0;
}
