#ifndef DUMPBBCVERTEXMAP_H__
#define DUMPBBCVERTEXMAP_H__

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

#endif /* __DUMPBBCVERTEXMAP_H__ */

