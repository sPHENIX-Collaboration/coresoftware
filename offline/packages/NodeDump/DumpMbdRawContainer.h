#ifndef NODEDUMP_DUMPMBDRAWCONTAINER_H
#define NODEDUMP_DUMPMBDRAWCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpMbdRawContainer : public DumpObject
{
 public:
  explicit DumpMbdRawContainer(const std::string &NodeName);
  virtual ~DumpMbdRawContainer() {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
