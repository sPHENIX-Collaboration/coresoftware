#ifndef NODEDUMP_DUMPBBCVERTEXMAP_H
#define NODEDUMP_DUMPBBCVERTEXMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpBbcVertexMap : public DumpObject
{
 public:
  DumpBbcVertexMap(const std::string &NodeName);
  virtual ~DumpBbcVertexMap() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif

