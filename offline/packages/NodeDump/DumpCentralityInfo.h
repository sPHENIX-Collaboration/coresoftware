#ifndef NODEDUMP_DUMPCENTRALITYINFO_H
#define NODEDUMP_DUMPCENTRALITYINFO_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpCentralityInfo : public DumpObject
{
 public:
  explicit DumpCentralityInfo(const std::string &NodeName);
  ~DumpCentralityInfo() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
