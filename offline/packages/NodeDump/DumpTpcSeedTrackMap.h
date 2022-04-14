#ifndef NODEDUMP_DUMPTPCSEEDTRACKMAP_H
#define NODEDUMP_DUMPTPCSEEDTRACKMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTpcSeedTrackMap : public DumpObject
{
 public:
  DumpTpcSeedTrackMap(const std::string &NodeName);
  ~DumpTpcSeedTrackMap() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
