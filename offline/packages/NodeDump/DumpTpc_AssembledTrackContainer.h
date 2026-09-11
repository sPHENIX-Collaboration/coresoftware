#ifndef NODEDUMP_DUMPTPCASSEMBLEDTRACKCONTAINER_H
#define NODEDUMP_DUMPTPCASSEMBLEDTRACKCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTpc_AssembledTrackContainer : public DumpObject
{
 public:
  explicit DumpTpc_AssembledTrackContainer(const std::string &NodeName);
  ~DumpTpc_AssembledTrackContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
