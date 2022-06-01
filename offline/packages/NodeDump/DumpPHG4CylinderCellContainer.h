#ifndef NODEDUMP_DUMPPHG4CYLINDERCELLCONTAINER_H
#define NODEDUMP_DUMPPHG4CYLINDERCELLCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4CylinderCellContainer : public DumpObject
{
 public:
  explicit DumpPHG4CylinderCellContainer(const std::string &NodeName);
  ~DumpPHG4CylinderCellContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
