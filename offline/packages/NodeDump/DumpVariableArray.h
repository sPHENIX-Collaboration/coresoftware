#ifndef NODEDUMP_DUMPVARIABLEARRAY_H
#define NODEDUMP_DUMPVARIABLEARRAY_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpVariableArray : public DumpObject
{
 public:
  explicit DumpVariableArray(const std::string &NodeName);
  ~DumpVariableArray() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
