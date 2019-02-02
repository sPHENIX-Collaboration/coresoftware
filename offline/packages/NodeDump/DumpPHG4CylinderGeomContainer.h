#ifndef NODEDUMP_DUMPPHG4CYLINDERGEOMCONTAINER_H
#define NODEDUMP_DUMPPHG4CYLINDERGEOMCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4CylinderGeomContainer : public DumpObject
{
 public:
  DumpPHG4CylinderGeomContainer(const std::string &NodeName);
  virtual ~DumpPHG4CylinderGeomContainer() {}

 protected:
  int process_Node(PHNode *mynode);
};

#endif
