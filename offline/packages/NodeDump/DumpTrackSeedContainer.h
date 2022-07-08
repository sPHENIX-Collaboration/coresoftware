#ifndef NODEDUMP_DUMPTRACKSEEDCONTAINER_H
#define NODEDUMP_DUMPTRACKSEEDCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTrackSeedContainer : public DumpObject
{
 public:
  explicit DumpTrackSeedContainer(const std::string &NodeName);
  ~DumpTrackSeedContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
