#ifndef NODEDUMP_DUMPMBDOUT_H
#define NODEDUMP_DUMPMBDOUT_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpMbdOut : public DumpObject
{
 public:
  explicit DumpMbdOut(const std::string &NodeName);
  virtual ~DumpMbdOut() {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
