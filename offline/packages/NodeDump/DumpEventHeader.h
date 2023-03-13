#ifndef NODEDUMP_DUMPEVENTHEADER_H
#define NODEDUMP_DUMPEVENTHEADER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpEventHeader : public DumpObject
{
 public:
  explicit DumpEventHeader(const std::string &NodeName);
  ~DumpEventHeader() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
