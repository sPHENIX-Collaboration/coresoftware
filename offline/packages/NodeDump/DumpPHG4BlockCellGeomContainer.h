#ifndef NODEDUMP_DUMPPHG4BLOCKCELLGEOMCONTAINER_H
#define NODEDUMP_DUMPPHG4BLOCKCELLGEOMCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4BlockCellGeomContainer : public DumpObject
{
 public:
  explicit DumpPHG4BlockCellGeomContainer(const std::string &NodeName);
  ~DumpPHG4BlockCellGeomContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
