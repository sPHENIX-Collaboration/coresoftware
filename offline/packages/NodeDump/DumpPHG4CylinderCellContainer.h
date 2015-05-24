#ifndef DUMPPHG4CYLINDERCELLCONTAINER_H__
#define DUMPPHG4CYLINDERCELLCONTAINER_H__

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4CylinderCellContainer : public DumpObject
{
 public:
  DumpPHG4CylinderCellContainer(const std::string &NodeName);
  virtual ~DumpPHG4CylinderCellContainer() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif /* DUMPPHG4CYLINDERCELLCONTAINER_H__ */

