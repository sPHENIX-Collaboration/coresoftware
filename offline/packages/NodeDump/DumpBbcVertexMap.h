#ifndef NODEDUMP_DUMPBBCVERTEXMAP_H
#define NODEDUMP_DUMPBBCVERTEXMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpBbcVertexMap : public DumpObject
{
 public:
  explicit DumpBbcVertexMap(const std::string &NodeName);
  ~DumpBbcVertexMap() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
