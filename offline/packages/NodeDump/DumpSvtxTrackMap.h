#ifndef NODEDUMP_DUMPSVTXTRACKMAP_H
#define NODEDUMP_DUMPSVTXTRACKMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpSvtxTrackMap : public DumpObject
{
 public:
  explicit DumpSvtxTrackMap(const std::string &NodeName);
  ~DumpSvtxTrackMap() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
