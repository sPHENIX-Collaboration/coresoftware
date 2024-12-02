#ifndef NODEDUMP_DUMPBBCPMTINFOCONTAINER_H
#define NODEDUMP_DUMPBBCPMTINFOCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpBbcPmtInfoContainer : public DumpObject
{
 public:
  explicit DumpBbcPmtInfoContainer(const std::string &NodeName);
  virtual ~DumpBbcPmtInfoContainer() {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
