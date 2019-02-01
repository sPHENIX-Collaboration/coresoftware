#ifndef NODEDUMP_DUMPSVTXCLUSTERMAP_H
#define NODEDUMP_DUMPSVTXCLUSTERMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpSvtxClusterMap : public DumpObject
{
 public:
  DumpSvtxClusterMap(const std::string &NodeName);
  virtual ~DumpSvtxClusterMap() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif

