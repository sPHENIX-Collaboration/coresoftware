#ifndef NODEDUMP_DUMPRAWTOWERCONTAINER_H
#define NODEDUMP_DUMPRAWTOWERCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpRawTowerContainer : public DumpObject
{
 public:
  DumpRawTowerContainer(const std::string &NodeName);
  virtual ~DumpRawTowerContainer() {}

 protected:
  int process_Node(PHNode *mynode);
};

#endif
