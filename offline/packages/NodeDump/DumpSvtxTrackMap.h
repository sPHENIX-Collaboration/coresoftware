#ifndef DUMPSVTXTRACKMAP_H__
#define DUMPSVTXTRACKMAP_H__

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpSvtxTrackMap : public DumpObject
{
 public:
  DumpSvtxTrackMap(const std::string &NodeName);
  virtual ~DumpSvtxTrackMap() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif /* DUMPSVTXTRACKMAP_H__ */

