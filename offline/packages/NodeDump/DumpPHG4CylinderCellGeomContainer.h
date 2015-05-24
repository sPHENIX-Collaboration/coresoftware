#ifndef DUMPPHG4CYLINDERCELLGEOMCONTAINER_H__
#define DUMPPHG4CYLINDERCELLGEOMCONTAINER_H__

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

#endif /* DUMPPHG4CYLINDERCELLGEOMCONTAINER_H__ */

