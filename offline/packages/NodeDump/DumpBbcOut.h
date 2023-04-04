#ifndef NODEDUMP_DUMPBBCOUT_H
#define NODEDUMP_DUMPBBCOUT_H


#include "DumpObject.h"

#include <string>

class PHNode;

class DumpBbcOut : public DumpObject
{
 public:
  DumpBbcOut(const std::string &NodeName);
  virtual ~DumpBbcOut() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif

