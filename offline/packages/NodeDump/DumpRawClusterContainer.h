#ifndef NODEDUMP_DUMPRAWCLUSTERCONTAINER_H
#define NODEDUMP_DUMPRAWCLUSTERCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpRawClusterContainer : public DumpObject
{
 public:
  explicit DumpRawClusterContainer(const std::string &NodeName);
  ~DumpRawClusterContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
