#ifndef NODEDUMP_DUMPRAWTOWERGEOMCONTAINER_H
#define NODEDUMP_DUMPRAWTOWERGEOMCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpRawTowerGeomContainer : public DumpObject
{
 public:
  explicit DumpRawTowerGeomContainer(const std::string &NodeName);
  ~DumpRawTowerGeomContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
