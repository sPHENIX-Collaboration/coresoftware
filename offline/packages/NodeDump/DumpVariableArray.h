#ifndef __DUMPVARIABLEARRAY_H__
#define __DUMPVARIABLEARRAY_H__

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

#endif /* __DUMPVARIABLEARRAY_H__ */

