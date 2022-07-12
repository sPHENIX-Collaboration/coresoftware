#ifndef NODEDUMP_DUMPSVTXVERTEXMAP_H
#define NODEDUMP_DUMPSVTXVERTEXMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpSvtxVertexMap : public DumpObject
{
 public:
  explicit DumpSvtxVertexMap(const std::string &NodeName);
  ~DumpSvtxVertexMap() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
