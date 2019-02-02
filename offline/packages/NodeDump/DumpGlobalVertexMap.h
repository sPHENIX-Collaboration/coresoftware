#ifndef NODEDUMP_DUMPGLOBALVERTEXMAP_H
#define NODEDUMP_DUMPGLOBALVERTEXMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpGlobalVertexMap : public DumpObject
{
 public:
  DumpGlobalVertexMap(const std::string &NodeName);
  virtual ~DumpGlobalVertexMap() {}

 protected:
  int process_Node(PHNode *mynode);
};

#endif
