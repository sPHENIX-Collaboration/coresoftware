#ifndef DUMPPHG4CYLINDERGEOMCONTAINER_H__
#define DUMPPHG4CYLINDERGEOMCONTAINER_H__

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

#endif /* DUMPPHG4CYLINDERGEOMCONTAINER_H__ */

