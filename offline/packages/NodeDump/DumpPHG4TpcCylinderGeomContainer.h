#ifndef NODEDUMP_DUMPPHG4TPCCYLINDERGEOMCONTAINER_H
#define NODEDUMP_DUMPPHG4TPCCYLINDERGEOMCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4TpcCylinderGeomContainer : public DumpObject
{
 public:
  explicit DumpPHG4TpcCylinderGeomContainer(const std::string &NodeName);
  ~DumpPHG4TpcCylinderGeomContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
