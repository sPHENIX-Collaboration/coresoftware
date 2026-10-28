#ifndef NODEDUMP_DUMPTRIGGERRUNINFO_H
#define NODEDUMP_DUMPTRIGGERRUNINFO_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTriggerRunInfo : public DumpObject
{
 public:
  explicit DumpTriggerRunInfo(const std::string &NodeName);
  ~DumpTriggerRunInfo() override {}

 protected:
  int process_Node(PHNode *mynode) override;
  int node_written = 0;
};

#endif
