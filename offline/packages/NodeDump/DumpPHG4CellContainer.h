#ifndef NODEDUMP_DUMPPHG4CELLCONTAINER_H
#define NODEDUMP_DUMPPHG4CELLCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4CellContainer : public DumpObject
{
 public:
  explicit DumpPHG4CellContainer(const std::string &NodeName);
  ~DumpPHG4CellContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
