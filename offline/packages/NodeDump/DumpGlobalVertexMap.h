#ifndef DUMPGLOBALVERTEXMAP_H__
#define DUMPGLOBALVERTEXMAP_H__

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

#endif /* __DUMPGLOBALVERTEXMAP_H__ */

