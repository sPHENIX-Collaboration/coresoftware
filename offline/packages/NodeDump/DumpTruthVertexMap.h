#ifndef NODEDUMP_DUMPTRUTHVERTEXMAP_H
#define NODEDUMP_DUMPTRUTHVERTEXMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTruthVertexMap : public DumpObject
{
 public:
  explicit DumpTruthVertexMap(const std::string &NodeName);
  ~DumpTruthVertexMap() override = default;

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
