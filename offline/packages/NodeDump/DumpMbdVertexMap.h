#ifndef NODEDUMP_DUMPMBDVERTEXMAP_H
#define NODEDUMP_DUMPMBDVERTEXMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpMbdVertexMap : public DumpObject
{
 public:
  explicit DumpMbdVertexMap(const std::string &NodeName);
  ~DumpMbdVertexMap() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
