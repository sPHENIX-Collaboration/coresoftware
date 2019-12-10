#ifndef NODEDUMP_DUMPSVTXTRACKMAP_H
#define NODEDUMP_DUMPSVTXTRACKMAP_H

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

#endif
