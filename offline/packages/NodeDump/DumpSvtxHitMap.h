#ifndef NODEDUMP_DUMPSVTXHITMAP_H
#define NODEDUMP_DUMPSVTXHITMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpSvtxHitMap : public DumpObject
{
 public:
  DumpSvtxHitMap(const std::string &NodeName);
  virtual ~DumpSvtxHitMap() {}

 protected:
  int process_Node(PHNode *mynode);
};

#endif
