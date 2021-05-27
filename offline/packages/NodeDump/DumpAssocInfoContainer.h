#ifndef NODEDUMP_DUMPASSOCINFOCONTAINER_H
#define NODEDUMP_DUMPASSOCINFOCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpAssocInfoContainer : public DumpObject
{
 public:
  DumpAssocInfoContainer(const std::string &NodeName);
  ~DumpAssocInfoContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
