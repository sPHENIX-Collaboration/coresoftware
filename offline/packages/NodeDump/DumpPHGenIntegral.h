#ifndef NODEDUMP_DUMPPHGENINTEGRAL_H
#define NODEDUMP_DUMPPHGENINTEGRAL_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHGenIntegral : public DumpObject
{
 public:
  explicit DumpPHGenIntegral(const std::string &NodeName);
  ~DumpPHGenIntegral() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
