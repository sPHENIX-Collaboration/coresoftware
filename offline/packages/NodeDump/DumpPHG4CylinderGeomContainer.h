#ifndef NODEDUMP_DUMPPHG4CYLINDERGEOMCONTAINER_H
#define NODEDUMP_DUMPPHG4CYLINDERGEOMCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4CylinderGeomContainer : public DumpObject
{
 public:
  explicit DumpPHG4CylinderGeomContainer(const std::string &NodeName);
  ~DumpPHG4CylinderGeomContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
