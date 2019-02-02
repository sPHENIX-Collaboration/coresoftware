#ifndef NODEDUMP_DUMPSVTXVERTEXMAP_H
#define NODEDUMP_DUMPSVTXVERTEXMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpSvtxVertexMap : public DumpObject
{
 public:
  DumpSvtxVertexMap(const std::string &NodeName);
  virtual ~DumpSvtxVertexMap() {}

 protected:
  int process_Node(PHNode *mynode);
};

#endif
