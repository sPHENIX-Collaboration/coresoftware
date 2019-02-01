#ifndef NODEDUMP_DUMPPHG4HITCONTAINER_H
#define NODEDUMP_DUMPPHG4HITCONTAINER_H

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

#endif

