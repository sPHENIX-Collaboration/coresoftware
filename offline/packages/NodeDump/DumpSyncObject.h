#ifndef NODEDUMP_DUMPSYNCOBJECT_H
#define NODEDUMP_DUMPSYNCOBJECT_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpSyncObject : public DumpObject
{
 public:
  explicit DumpSyncObject(const std::string &NodeName);
  ~DumpSyncObject() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
