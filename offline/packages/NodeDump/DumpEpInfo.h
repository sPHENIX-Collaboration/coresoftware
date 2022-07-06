#ifndef NODEDUMP_DUMPEPINFO_H
#define NODEDUMP_DUMPEPINFO_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpEpInfo : public DumpObject
{
 public:
  explicit DumpEpInfo(const std::string &NodeName);
  ~DumpEpInfo() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
