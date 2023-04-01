#ifndef NODEDUMP_DUMPBBCPMTCONTAINER_H
#define NODEDUMP_DUMPBBCPMTCONTAINER_H


#include "DumpObject.h"

#include <string>

class PHNode;

class DumpBbcPmtContainer : public DumpObject
{
 public:
  DumpBbcPmtContainer(const std::string &NodeName);
  virtual ~DumpBbcPmtContainer() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif

