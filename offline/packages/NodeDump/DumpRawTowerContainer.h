#ifndef NODEDUMP_DUMPRAWTOWERCONTAINER_H
#define NODEDUMP_DUMPRAWTOWERCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpRawTowerContainer : public DumpObject
{
 public:
  DumpRawTowerContainer(const std::string &NodeName);
  ~DumpRawTowerContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
