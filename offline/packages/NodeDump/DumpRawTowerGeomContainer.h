#ifndef NODEDUMP_DUMPRAWTOWERGEOMCONTAINER_H
#define NODEDUMP_DUMPRAWTOWERGEOMCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpRawTowerGeomContainer : public DumpObject
{
 public:
  DumpRawTowerGeomContainer(const std::string &NodeName);
  virtual ~DumpRawTowerGeomContainer() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif

