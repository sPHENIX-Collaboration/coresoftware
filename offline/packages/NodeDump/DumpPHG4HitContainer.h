#ifndef DUMPPHG4HITCONTAINER_H__
#define DUMPPHG4HITCONTAINER_H__

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4HitContainer : public DumpObject
{
 public:
  DumpPHG4HitContainer(const std::string &NodeName);
  virtual ~DumpPHG4HitContainer() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif /* __DUMPPHG4HITCONTAINER_H__ */

