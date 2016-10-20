#ifndef DUMPSVTXHITMAP_H__
#define DUMPSVTXHITMAP_H__

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

#endif /* DUMPSVTXHITMAP_H__ */

