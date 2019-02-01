#ifndef NODEDUMP_DUMPRAWCLUSTERCONTAINER_H
#define NODEDUMP_DUMPRAWCLUSTERCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpRawClusterContainer : public DumpObject
{
 public:
  DumpRawClusterContainer(const std::string &NodeName);
  virtual ~DumpRawClusterContainer() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif

