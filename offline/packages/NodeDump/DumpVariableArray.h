#ifndef NODEDUMP_DUMPVARIABLEARRAY_H
#define NODEDUMP_DUMPVARIABLEARRAY_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpVariableArray : public DumpObject
{
 public:
  DumpVariableArray(const std::string &NodeName);
  virtual ~DumpVariableArray() {}

 protected:
  int process_Node(PHNode *mynode);
};

#endif
