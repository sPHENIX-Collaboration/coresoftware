#ifndef NODEDUMP_DUMPRAWTOWERGEOM_H
#define NODEDUMP_DUMPRAWTOWERGEOM_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpRawTowerGeom : public DumpObject
{
 public:
  DumpRawTowerGeom(const std::string &NodeName);
  ~DumpRawTowerGeom() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
