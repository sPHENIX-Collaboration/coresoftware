#ifndef NODEDUMP_DUMPPHG4BLOCKGEOMCONTAINER_H
#define NODEDUMP_DUMPPHG4BLOCKGEOMCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4BlockGeomContainer : public DumpObject
{
 public:
  explicit DumpPHG4BlockGeomContainer(const std::string &NodeName);
  ~DumpPHG4BlockGeomContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
