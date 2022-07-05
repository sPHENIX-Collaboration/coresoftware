#ifndef NODEDUMP_DUMPGLOBALVERTEXMAP_H
#define NODEDUMP_DUMPGLOBALVERTEXMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpGlobalVertexMap : public DumpObject
{
 public:
  explicit DumpGlobalVertexMap(const std::string &NodeName);
  ~DumpGlobalVertexMap() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
