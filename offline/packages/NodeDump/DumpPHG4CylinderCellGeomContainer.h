#ifndef NODEDUMP_DUMPPHG4CYLINDERCELLGEOMCONTAINER_H
#define NODEDUMP_DUMPPHG4CYLINDERCELLGEOMCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4CylinderCellGeomContainer : public DumpObject
{
 public:
  DumpPHG4CylinderCellGeomContainer(const std::string &NodeName);
  virtual ~DumpPHG4CylinderCellGeomContainer() {}

 protected:
  int process_Node(PHNode *mynode);
};

#endif
