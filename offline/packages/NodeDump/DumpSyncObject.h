#ifndef NODEDUMP_DUMPSYNCOBJECT_H
#define NODEDUMP_DUMPSYNCOBJECT_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpSyncObject : public DumpObject
{
 public:
  DumpSyncObject(const std::string &NodeName);
  virtual ~DumpSyncObject() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif

