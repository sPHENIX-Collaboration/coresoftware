#ifndef DUMPSVTXVERTEXMAP_H__
#define DUMPSVTXVERTEXMAP_H__

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

#endif /* DUMPSVTXVERTEXMAP_H__ */

