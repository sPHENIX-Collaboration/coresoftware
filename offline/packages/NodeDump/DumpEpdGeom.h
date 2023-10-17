#ifndef NODEDUMP_DUMPEPDGEOM_H
#define NODEDUMP_DUMPEPDGEOM_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpEpdGeom : public DumpObject
{
 public:
  explicit DumpEpdGeom(const std::string &NodeName);
  ~DumpEpdGeom() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
