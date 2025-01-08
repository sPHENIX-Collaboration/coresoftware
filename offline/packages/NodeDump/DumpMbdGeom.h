#ifndef NODEDUMP_DUMPMBDGEOM_H
#define NODEDUMP_DUMPMBDGEOM_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpMbdGeom : public DumpObject
{
 public:
  explicit DumpMbdGeom(const std::string &NodeName);
  ~DumpMbdGeom() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
