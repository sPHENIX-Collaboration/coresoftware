#ifndef NODEDUMP_DUMPMBDPMTCONTAINER_H
#define NODEDUMP_DUMPMBDPMTCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpMbdPmtContainer : public DumpObject
{
 public:
  explicit DumpMbdPmtContainer(const std::string &NodeName);
  virtual ~DumpMbdPmtContainer() {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
