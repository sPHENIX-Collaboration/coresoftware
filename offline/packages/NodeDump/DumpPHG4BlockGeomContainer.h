#ifndef NODEDUMP_DUMPPHG4BLOCKGEOMCONTAINER_H
#define NODEDUMP_DUMPPHG4BLOCKGEOMCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4BlockGeomContainer : public DumpObject
{
 public:
  DumpPHG4BlockGeomContainer(const std::string &NodeName);
  virtual ~DumpPHG4BlockGeomContainer() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif

