#ifndef __DUMPSYNCOBJECT_H__
#define __DUMPSYNCOBJECT_H__

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

#endif /* __DUMPSYNCOBJECT_H__ */

