#ifndef NODEDUMP_DUMPFLAGSAVE_H
#define NODEDUMP_DUMPFLAGSAVE_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpFlagSave : public DumpObject
{
 public:
  explicit DumpFlagSave(const std::string &NodeName);
  ~DumpFlagSave() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
