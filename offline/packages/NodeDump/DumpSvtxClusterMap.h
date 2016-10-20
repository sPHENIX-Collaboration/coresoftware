#ifndef DUMPSVTXCLUSTERMAP_H__
#define DUMPSVTXCLUSTERMAP_H__

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

#endif /* DUMPSVTXCLUSTERMAP_H__ */

