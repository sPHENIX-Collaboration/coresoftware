#ifndef NODEDUMP_DUMPPHG4CYLINDERCELLGEOMCONTAINER_H
#define NODEDUMP_DUMPPHG4CYLINDERCELLGEOMCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4CylinderCellGeomContainer : public DumpObject
{
 public:
  explicit DumpPHG4CylinderCellGeomContainer(const std::string &NodeName);
  ~DumpPHG4CylinderCellGeomContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
